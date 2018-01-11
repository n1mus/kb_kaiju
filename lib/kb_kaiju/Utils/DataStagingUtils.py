import os
import time
import glob
import subprocess

from Workspace.WorkspaceClient import Workspace
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIServiceClient import SetAPI


class DataStagingUtils(object):

    def __init__(self, config, ctx):
        self.ctx = ctx
        self.scratch = os.path.abspath(config['scratch'])
        self.ws_url = config['workspace-url']
        self.serviceWizardURL = config['srv-wiz-url']
        self.callbackURL = config['SDK_CALLBACK_URL']
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        self.SE_flag = 'SE'
        self.PE_lfag = 'PE'


    def stage_input(self, input_refs, fasta_file_extension):
        '''
        Stage input based on an input data reference for Kaiju

        input_refs can be a list of references to a PairedEndLibrary, a SingleEndLibrary, or a ReadsSet

        This method creates a directory in the scratch area with the set of Fasta/Fastq files, names
        will have the fasta_file_extension parameter tacked on.

            ex:

            staged_input = stage_input('124/15/1', 'fastq')

            staged_input
            {"input_dir": '...'}
        '''
        # config
        #SERVICE_VER = 'dev'
        SERVICE_VER = 'release'

        # generate a folder in scratch to hold the input
        suffix = str(int(time.time() * 1000))
        input_dir = os.path.join(self.scratch, 'input_reads_' + suffix)
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)

        # ruClient
        try:
            ruClient = ReadsUtils(self.callbackURL, token=self.ctx['token'], service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate ruClient with callbackURL: '+ self.callbackURL +' ERROR: ' + str(e))

        # setAPI_Client
        try:
            #setAPI_Client = SetAPI (url=self.callbackURL, token=self.ctx['token'])  # for SDK local.  local doesn't work for SetAPI
            setAPI_Client = SetAPI (url=self.serviceWizardURL, token=self.ctx['token'])  # for dynamic service
        except Exception as e:
            raise ValueError('Unable to instantiate setAPI_Client with serviceWizardURL: '+ self.serviceWizardURL +' ERROR: ' + str(e))


        # 2) expand any sets and build a non-redundant list of reads input objs
        ws = Workspace(self.ws_url)
        expanded_input = []
        input_ref_seen = dict()
        SE_types = ['KBaseFile.SingleEndLibrary', 'KBaseAssembly.SingleEndLibrary']
        PE_types = ['KBaseFile.PairedEndLibrary', 'KBaseAssembly.PairedEndLibrary']

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        for input_ref in input_refs:
            input_info = ws.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0]
            obj_name = input_info[NAME_I]
            type_name = input_info[TYPE_I].split('-')[0]

            # ReadsSet
            if type_name in ['KBaseSets.ReadsSet']:
                try:
                    input_readsSet_obj = setAPI_Client.get_reads_set_v1 ({'ref':input_ref,'include_item_info':1})

                except Exception as e:
                    raise ValueError('SetAPI FAILURE: Unable to get read library set object from workspace: (' + str(input_ref)+")\n" + str(e))

                for readsLibrary_obj in input_readsSet_obj['data']['items']:
                    this_reads_ref = readsLibrary_obj['ref']
                    if this_reads_ref in input_ref_seen:
                        continue
                    input_ref_seen[this_reads_ref] = True

                    this_reads_name = readsLibrary_obj['info'][NAME_I]
                    reads_item_type = readsLibrary_obj['info'][TYPE_I]
                    reads_item_type = re.sub ('-[0-9]+\.[0-9]+$', "", reads_item_type)  # remove trailing version
                    if reads_item_type in PE_types:
                        this_read_type = self.PE_flag
                    elif reads_item_type in SE_types:
                        this_read_type = self.SE_flag
                    else:
                        raise ValueError ("Can't handle read item type '"+reads_item_type+"' obj_name: '"+this_reads_name+" in Set: '"+str(input_ref)+"'")
                    expanded_input.append({'ref':  this_reads_ref,
                                           'name': this_reads_name,
                                           'type': this_reads_type
                                       })
            # SingleEnd Library
            elif type_name in [SE_types]:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.SE_flag
                expanded_input.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                   })
            # PairedEnd Library
            elif type_name in [PE_types]:
                this_reads_ref = input_ref
                if this_reads_ref in input_ref_seen:
                    continue
                input_ref_seen[this_reads_ref] = True
                this_reads_name = obj_name
                this_reads_type = self.PE_flag
                expanded_input.append({'ref':  this_reads_ref,
                                       'name': this_reads_name,
                                       'type': this_reads_type
                                   })
            else:
                raise ValueError ("Illegal type in input_refs: "+str(obj_name)+" ("+str(input_ref)+") is of type: '"+str(type_name)+"'")


        # 3) Download reads
        for input_item_i,input_item in enumerate(expanded_input):

            try:
                readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [input_item['ref']],
                                                                  'interleaved': 'false'})
            except Exception as e:
                raise ValueError('Unable to get read library object from workspace: (' + str(input_item['ref']) +")\n" + str(e))

            # PE Lib
            if input_item['type'] == self.PE_flag:
                input_fwd_file_path = readsLibrary['files'][input_item['ref']]['files']['fwd']
                input_rev_file_path = readsLibrary['files'][input_item['ref']]['files']['rev']
                fwd_filename = os.path.join(input_dir, input_item['name'] + '.fwd.' + fasta_file_extension)
                rev_filename = os.path.join(input_dir, input_item['name'] + '.rev.' + fasta_file_extension)
                if input_fwd_file_path != fwd_filename:
                    shutil.move(input_fwd_file_path, fwd_filename)
                if input_rev_file_path != rev_filename:
                    shutil.move(input_rev_file_path, rev_filename)
                expanded_input[input_item_i]['fwd_file'] = fwd_filename
                expanded_input[input_item_i]['rev_file'] = rev_filename

                if not os.path.isfile(fwd_filename):
                    raise ValueError('Error generating reads file '+fwd_filename)
                if not os.path.isfile(rev_filename):
                    raise ValueError('Error generating reads file '+rev_filename)
                # make sure fasta file isn't empty
                min_fasta_len = 1
                if not self.fasta_seq_len_at_least(fwd_filename, min_fasta_len):
                    raise ValueError('Reads Library is empty in filename: '+str(fwd_filename))
                if not self.fasta_seq_len_at_least(input_rev_file_path, min_fasta_len):
                    raise ValueError('Reads Library is empty in filename: '+str(rev_filename))

            elif input_item['type'] == self.SE_flag:
                input_fwd_file_path = readsLibrary['files'][input_item['ref']]['files']['fwd']
                fwd_filename = os.path.join(input_dir, input_item['name'] + '.fwd.' + fasta_file_extension)
                if input_fwd_file_path != fwd_filename:
                    shutil.move(input_fwd_file_path, fwd_filename)
                expanded_input[input_item_i]['fwd_file'] = fwd_filename

                if not os.path.isfile(fwd_filename):
                    raise ValueError('Error generating reads file '+fwd_filename)
                # make sure fasta file isn't empty
                min_fasta_len = 1
                if not self.fasta_seq_len_at_least(fwd_filename, min_fasta_len):
                    raise ValueError('Reads Library is empty in filename: '+str(fwd_filename))
            
            
        return {'input_dir': input_dir, 'folder_suffix': suffix, 'expanded_input': expanded_input}


    def fasta_seq_len_at_least(self, fasta_path, min_fasta_len=1):
        '''
        counts the number of non-header, non-whitespace characters in a FASTA file
        '''
        seq_len = 0
        with open (fasta_path, 'r', 0) as fasta_handle:
            for line in fasta_handle:
                line = line.strip()
                if line.startswith('>'):
                    continue
                line = line.replace(' ','')
                seq_len += len(line)
                if seq_len >= min_fasta_len:
                    return True
        return False
