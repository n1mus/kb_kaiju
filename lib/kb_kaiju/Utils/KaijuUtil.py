# -*- coding: utf-8 -*-
import time
import os
import uuid
import subprocess
import sys

from KBaseReport.KBaseReportClient import KBaseReport

from kb_kaiju.Utils.DataStagingUtils import DataStagingUtils
from kb_kaiju.Utils.OutputBuilder import OutputBuilder


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))
    sys.stdout.flush()


class KaijuUtil:

    def __init__(self, config, ctx):
        self.config = config
        self.ctx = ctx
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.threads = config['threads']
        self.SE_flag = 'SE'
        self.PE_flag = 'PE'

        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        

    def run_kaiju_with_krona(self, params):
        '''
        Main entry point for running kaiju + krona as a KBase App
        '''

        # 0) validate basic parameters
        method = 'run_kaiju_and_krona'
        required_params = ['workspace_name',
                           'input_refs',
                           'output_biom_name',
                           'tax_levels',
                           'db_type',
                           'seg_filter',
                           'min_match_length',
                           'greedy_run_mode',
                           'filter_percent',
                           'filter_unclassified',
                           'full_tax_path'
                          ]
        for arg in required_params:
            if arg not in params or params[arg] == None or params[arg] == '':
                raise ValueError ("Must define required param: '"+arg+"' for method: '"+str(method)+"()'")

        if 'greedy_run_mode' in params and int(params['greedy_run_mode']) == 1:
            greedy_required_params = ['greedy_allowed_mismatches',
                                      'greedy_min_match_score',
                                      'greedy_max_e_value'
                                  ]
            for arg in greedy_required_params:
                if arg not in params or params[arg] == None or params[arg] == '':
                    raise ValueError ("Must define GREEDY MODE required param: '"+arg+"' for method: '"+str(method)+"()'")


        # 1) stage input data
        dsu = DataStagingUtils(self.config, self.ctx)
        staged_input = dsu.stage_input(params['input_refs'], 'fastq')
        input_dir = staged_input['input_dir']
        suffix = staged_input['folder_suffix']
        expanded_input = staged_input['expanded_input']

        #print ("INPUT_DIR: "+input_dir)
        #print ("SUFFIX: "+suffix)
        #print ("EXPANDED_INPUT: ")
        #print (expanded_input)
            

        output_dir = os.path.join(self.scratch, 'output_' + suffix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        #plots_dir = os.path.join(self.scratch, 'plot_' + suffix)
        html_dir = os.path.join(self.scratch, 'html_' + suffix)
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)

        log('Staged input directory: ' + input_dir)


        # 2) run Kaiju in batch
        kaiju_output_folder = os.path.join(output_dir, 'kaiju_output')
        if not os.path.exists(kaiju_output_folder):
            os.makedirs(kaiju_output_folder)
        kaiju_options = {'input_reads':               expanded_input,
                         'out_folder':                kaiju_output_folder,
                         'tax_levels':                params['tax_levels'],
                         'db_type':                   params['db_type'],
                         'seg_filter':                params['seg_filter'],
                         'min_match_length':          params['min_match_length'],
                         'greedy_run_mode':           params['greedy_run_mode'],
                         'greedy_allowed_mismatches': params['greedy_allowed_mismatches'],
                         'greedy_min_match_score':    params['greedy_min_match_score'],
                         'threads':                   self.threads
                        }
        self.run_kaiju_batch (kaiju_options)


        # 3) create Summary Reports in batch
        kaijuReport_output_folder = os.path.join(output_dir, 'kaiju_report')
        if not os.path.exists(kaijuReport_output_folder):
            os.makedirs(kaijuReport_output_folder)
        kaijuReport_options = {'input_reads':               expanded_input,
                               'in_folder':                 kaiju_output_folder,
                               'out_folder':                kaijuReport_output_folder,
                               'tax_levels':                params['tax_levels'],
                               'db_type':                   params['db_type'],
                               'filter_percent':            params['filter_percent'],
                               'filter_unclassified':       params['filter_unclassified'],
                               'full_tax_path':             params['full_tax_path']
                           }
        self.run_kaijuReport_batch (kaijuReport_options)


        # 4) create HTML Summary Reports in batch
        kaijuReportHTML_options = {'input_reads':               expanded_input,
                                   'in_folder':                 kaijuReport_output_folder,
                                   'out_folder':                html_dir,
                                   'tax_levels':                params['tax_levels']
        }
        #self.run_kaijuReportHTML_batch (kaijuReportHTML_options)


        # 5) create Krona plots
        #self.build_checkM_lineage_wf_plots(input_dir, output_dir, plots_dir, all_seq_fasta_file, tetra_file)


        # 6) Package results
        #outputBuilder = OutputBuilder(output_dir, plots_dir, self.scratch, self.callback_url)
        #output_packages = self._build_output_packages(params, outputBuilder, input_dir)


        # 7) build the HTML report
        #os.makedirs(html_dir)
        #outputBuilder.build_html_output_for_lineage_wf(html_dir, params['input_ref'])
        #html_zipped = outputBuilder.package_folder(html_dir, 'report.html', 'Summarized report from CheckM')


        # 8) save report
        report_params = {'message': '',
                         #'direct_html_link_index': 0,
                         #'html_links': [html_zipped],
                         #'file_links': output_packages,
                         'report_object_name': 'kb_kaiju_report_' + str(uuid.uuid4()),
                         'workspace_name': params['workspace_name']
                         }

        kr = KBaseReport(self.callback_url)
        report_output = kr.create_extended_report(report_params)

        returnVal = {'report_name': report_output['name'],
                     'report_ref':  report_output['ref']}
        return returnVal


    def run_proc(self, command, log_output_file=None):
        log('Running: ' + ' '.join(command))
        
        if log_output_file:  # if output is too chatty for STDOUT
            open (log_output_file, 'w')
            p = subprocess.Popen(command, cwd=self.scratch, shell=False, stdout=log_output_file, stderr=subprocess.STDOUT)
        else:
            p = subprocess.Popen(command, cwd=self.scratch, shell=False)
        exitCode = p.wait()
            
        if log_output_file:
            log_output_file.close()

        if (exitCode == 0):
            log('Executed command: ' + ' '.join(command) + '\n' +
                'Exit Code: ' + str(exitCode))
        else:
            raise ValueError('Error running command: ' + ' '.join(command) + '\n' +
                             'Exit Code: ' + str(exitCode))
        return exitCode
            

    def run_kaiju_batch(self, options, dropOutput=False):
        input_reads = options['input_reads']
        for input_reads_item in input_reads:
            single_kaiju_run_options = options
            single_kaiju_run_options['input_item'] = input_reads_item

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaiju' + '.stdout')
            
            command = self._build_kaiju_command(single_kaiju_run_options)
            self.run_proc (command, log_output_file)


    def run_kaijuReport_batch(self, options, dropOutput=False):
        input_reads = options['input_reads']
        for input_reads_item in input_reads:
            for tax_level in options['tax_levels']:
                single_kaijuReport_run_options = options
                single_kaijuReport_run_options['input_item'] = input_reads_item
                single_kaijuReport_run_options['tax_level'] = tax_level
            
                log_output_file = None
                if dropOutput:  # if output is too chatty for STDOUT
                    log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaijuReport' + '.stdout')
            
                command = self._build_kaijuReport_command(single_kaijuReport_run_options)
                self.run_proc (command, log_output_file)


    def _validate_kaiju_options(self, options):
        # 1st order required
        func_name = 'kaiju'
        required_opts = [ 'db_type',
                          'min_match_length',
                          'greedy_run_mode'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # 2nd order required
        if 'greedy_run_mode' in options and int(options['greedy_run_mode']) == 1:
            opt = 'greedy_allowed_mismatches'
            if opt not in options or int(options[opt]) < 1:
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()' if running in greedy_run_mode")

        # input file validation
        if not os.path.getsize(options['input_item']['fwd_file']) > 0:
            raise ValueError ('missing or empty fwd reads file: '+options['input_item']['fwd_file'])
        if options['input_item']['type'] == self.PE_flag:
            if not os.path.getsize(options['input_item']['rev_file']) > 0:
                raise ValueError ('missing or empty rev reads file: '+options['input_item']['rev_file'])

        # db validation
        DB = 'KAIJU_DB_PATH'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaiju_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_PATH'):
            command_list.append('-f')
            command_list.append(str(options.get('KAIJU_DB_PATH')))
        if options['input_item'].get('fwd_file'):
            command_list.append('-i')
            command_list.append(str(options['input_item'].get('fwd_file')))
        if options['input_item'].get('type') == self.PE_flag:
            command_list.append('-j')
            command_list.append(str(options['input_item'].get('rev_file')))
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'.kaiju'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if int(options.get('seg_filter')) == 1:
            command_list.append('-x')
        if options.get('min_match_length'):
            command_list.append('-m')
            command_list.append(str(options.get('min_match_length')))
        if int(options.get('greedy_run_mode')) == 1:
            command_list.append('-a')
            command_list.append('greedy')
            if options.get('greedy_allowed_mismatches'):
                command_list.append('-e')
                command_list.append(str(options.get('greedy_allowed_mismatches')))
            if options.get('greedy_min_match_score'):
                command_list.append('-s')
                command_list.append(str(options.get('greedy_min_match_score')))

        if options.get('threads'):
            command_list.append('-z')
            command_list.append(str(options.get('threads')))
        if options.get('verbose'):
            command_list.append('-v')
        

    def _build_kaiju_command(self, options, verbose=True):

        KAIJU_BIN_DIR  = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_BIN      = os.path.join(KAIJU_BIN_DIR, 'kaiju')
        KAIJU_DB_DIR   = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['verbose'] = verbose
        if self.threads and self.threads > 1:
            options['threads'] = self.threads

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        #options['KAIJU_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')  # don't need for kaiju cmd

        if options['db_type'] == 'kaiju_index':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_pg':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_nr':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr.fmi')
        elif options['db_type'] == 'kaiju_index_nr_euk':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr_euk.fmi')
        else:
            raise ValueError ('bad db_type: '+options['db_type']+' (must be one of "kaiju_index", "kaiju_index_pg", "kaiju_index_nr", "kaiju_index_nr_euk")')

        self._validate_kaiju_options(options)

        command = [KAIJU_BIN]
        self._process_kaiju_options(command, options)

        return command



    def _validate_kaijuReport_options(self, options):
        # 1st order required
        func_name = 'kaijuReport'
        required_opts = [ 'db_type',

                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # 2nd order required
        if 'greedy_run_mode' in options and int(options['greedy_run_mode']) == 1:
            opt = 'greedy_allowed_mismatches'
            if opt not in options or int(options[opt]) < 1:
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()' if running in greedy_run_mode")

        # input file validation
        if not os.path.getsize(options['input_item']['fwd_file']) > 0:
            raise ValueError ('missing or empty fwd reads file: '+options['input_item']['fwd_file'])
        if options['input_item']['type'] == self.PE_flag:
            if not os.path.getsize(options['input_item']['rev_file']) > 0:
                raise ValueError ('missing or empty rev reads file: '+options['input_item']['rev_file'])

        # db validation
        DB = 'KAIJU_DB_PATH'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaijuReport_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_PATH'):
            command_list.append('-f')
            command_list.append(str(options.get('KAIJU_DB_PATH')))
        if options['input_item'].get('fwd_file'):
            command_list.append('-i')
            command_list.append(str(options['input_item'].get('fwd_file')))
        if options['input_item'].get('type') == self.PE_flag:
            command_list.append('-j')
            command_list.append(str(options['input_item'].get('rev_file')))
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'.kaiju'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if int(options.get('seg_filter')) == 1:
            command_list.append('-x')
        if options.get('min_match_length'):
            command_list.append('-m')
            command_list.append(str(options.get('min_match_length')))
        if int(options.get('greedy_run_mode')) == 1:
            command_list.append('-a')
            command_list.append('greedy')
            if options.get('greedy_allowed_mismatches'):
                command_list.append('-e')
                command_list.append(str(options.get('greedy_allowed_mismatches')))
            if options.get('greedy_min_match_score'):
                command_list.append('-s')
                command_list.append(str(options.get('greedy_min_match_score')))

        if options.get('threads'):
            command_list.append('-z')
            command_list.append(str(options.get('threads')))
        if options.get('verbose'):
            command_list.append('-v')
        

    def _build_kaijuReport_command(self, options, verbose=True):

        KAIJU_BIN_DIR  = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_BIN      = os.path.join(KAIJU_BIN_DIR, 'kaiju')
        KAIJU_DB_DIR   = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['verbose'] = verbose
        if self.threads and self.threads > 1:
            options['threads'] = self.threads

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        #options['KAIJU_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')  # don't need for kaiju cmd

        if options['db_type'] == 'kaiju_index':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_pg':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db.fmi')
        elif options['db_type'] == 'kaiju_index_nr':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr.fmi')
        elif options['db_type'] == 'kaiju_index_nr_euk':
            options['KAIJU_DB_PATH'] = os.path.join(KAIJU_DB_DIR, 'kaiju_db_nr_euk.fmi')
        else:
            raise ValueError ('bad db_type: '+options['db_type']+' (must be one of "kaiju_index", "kaiju_index_pg", "kaiju_index_nr", "kaiju_index_nr_euk")')

        self._validate_kaiju_options(options)

        command = [KAIJU_BIN]
        self._process_kaiju_options(command, options)

        return command


    def _build_output_packages(self, params, outputBuilder, input_dir):

        output_packages = []

        #if 'save_output_dir' in params and str(params['save_output_dir']) == '1':
        if True:
            log('packaging full output directory')
            zipped_output_file = outputBuilder.package_folder(outputBuilder.output_dir, 'full_output.zip',
                                                              'Full output of CheckM')
            output_packages.append(zipped_output_file)
        else:  # ADD LATER?
            log('not packaging full output directory, selecting specific files')
            crit_out_dir = os.path.join(self.scratch, 'critical_output_' + os.path.basename(input_dir))
            os.makedirs(crit_out_dir)
            zipped_output_file = outputBuilder.package_folder(outputBuilder.output_dir, 'selected_output.zip',
                                                              'Selected output from the CheckM analysis')
            output_packages.append(zipped_output_file)


        if 'save_plots_dir' in params and str(params['save_plots_dir']) == '1':
            log('packaging output plots directory')
            zipped_output_file = outputBuilder.package_folder(outputBuilder.plots_dir, 'plots.zip',
                                                              'Output plots from CheckM')
            output_packages.append(zipped_output_file)
        else:
            log('not packaging output plots directory')

        return output_packages
