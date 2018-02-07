import os
import shutil
import ast
import sys
import time
import re

import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random
from random import shuffle

from biokbase.workspace.client import Workspace as workspaceService
#from Workspace.WorkspaceClient import Workspace as workspaceService
from DataFileUtil.DataFileUtilClient import DataFileUtil


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))
    sys.stdout.flush()


class OutputBuilder(object):
    '''
    Constructs the output HTML report and artifacts based on Kaiju and Krona
    runs.  This includes creating matplotlib graphs of the summaries and
    modifying the Krona HTML to offer tabbed href links between html pages
    '''

    def __init__(self, output_folders, scratch_dir, callback_url, workspace_url):
        self.output_folders = output_folders
        self.scratch = scratch_dir
        self.callback_url = callback_url
        self.workspace_url = workspace_url
        self.wsClient = None

        # store parsed info
        self.parsed_summary = dict()

        # leave out light colors
        self.no_light_color_names = [
            #'aliceblue',
            'aqua',
            'aquamarine',
            #'azure',
            #'beige',
            #'bisque',
            #'blanchedalmond',
            'blue',
            'blueviolet',
            'brown',
            'burlywood',
            'cadetblue',
            'chartreuse',
            'chocolate',
            'coral',
            'cornflowerblue',
            #'cornsilk',
            'crimson',
            'cyan',
            'darkblue',
            'darkcyan',
            'darkgoldenrod',
            'darkgreen',
            'darkkhaki',
            'darkmagenta',
            'darkolivegreen',
            'darkorange',
            'darkorchid',
            'darkred',
            'darksalmon',
            'darkseagreen',
            'darkslateblue',
            #'darkslategray',
            'darkturquoise',
            'darkviolet',
            'deeppink',
            'deepskyblue',
            'dodgerblue',
            'firebrick',
            'forestgreen',
            'fuchsia',
            #'gainsboro',
            'gold',
            'goldenrod',
            'green',
            'greenyellow',
            #'honeydew',
            'hotpink',
            'indianred',
            'indigo',
            'khaki',
            #'lavender',
            #'lavenderblush',
            'lawngreen',
            #'lemonchiffon',
            'lightblue',
            #'lightcoral',
            #'lightcyan'
            #'lightgoldenrodyellow',
            'lightgreen',
            'lightpink',
            'lightsalmon',
            'lightseagreen',
            'lightskyblue',
            #'lightslategray',
            #'lightsteelblue',
            #'lightyellow',
            'lime',
            'limegreen',
            #'magenta',   # magenta reserved for viruses
            'maroon',
            'mediumaquamarine',
            'mediumblue',
            'mediumorchid',
            'mediumpurple',
            'mediumseagreen',
            'mediumslateblue',
            'mediumspringgreen',
            'mediumturquoise',
            'mediumvioletred',
            'midnightblue',
            #'mintcream',
            #'mistyrose',
            #'moccasin',
            'navy',
            #'oldlace',
            'olive',
            'olivedrab',
            'orange',
            'orangered',
            'orchid',
            #'palegoldenrod',
            'palegreen',
            'paleturquoise',
            'palevioletred',
            #'papayawhip',
            'peachpuff',
            #'peru',
            'pink',
            'plum',
            'powderblue',
            'purple',
            'red',
            'rosybrown',
            'royalblue',
            'saddlebrown',
            'salmon',
            'sandybrown',
            'seagreen',
            #'seashell',
            'sienna',
            'skyblue',
            'slateblue',
            'springgreen',
            'steelblue',
            #'tan',
            'teal',
            #'thistle',
            'tomato',
            'turquoise',
            'violet',
            #'wheat',
            #'yellow',
            #'yellowgreen'
        ]


    def generate_sparse_biom1_0_matrix(self, ctx, options): 
        tax_level       = options['tax_level']
        input_reads     = options['input_reads']
        in_folder       = options['in_folder']
        workspace_name  = options['workspace_name']
        output_obj_name = options['output_obj_name']

        abundance_matrix = []
        abundance_by_sample = []
        lineage_seen = dict()
        lineage_order = []
        #extra_bucket_order = []
        sample_order = []
        #classified_frac = []
        biom_obj = dict()

        # parse summary
        for input_reads_item in input_reads:
            sample_order.append(input_reads_item['name'])

            this_summary_file = os.path.join (in_folder, input_reads_item['name']+'-'+tax_level+'.kaijuReport')
            (this_abundance, this_lineage_order, this_classified_frac) = self._parse_kaiju_summary_file (this_summary_file, tax_level)
            for lineage_name in this_lineage_order:
                if lineage_name not in lineage_seen:
                    lineage_seen[lineage_name] = True
                    if lineage_name.startswith('tail (<') \
                       or lineage_name.startswith('viruses') \
                       or lineage_name.startswith('unassigned at'):
                        #extra_bucket_order.append(lineage_name)
                        continue
                    else:
                        lineage_order.append(lineage_name)
            abundance_by_sample.append(this_abundance)
            #classified_frac.append(this_classified_frac)

        # create sparse matrix (note: vals in each sample do not sum to 100% because we're dumping buckets)
        biom_data = []
        for lineage_i,lineage_name in enumerate(lineage_order):
            for sample_i,sample_name in enumerate(sample_order):
                if lineage_name in abundance_by_sample[sample_i]:
                    biom_data.append([lineage_i, sample_i, abundance_by_sample[sample_i][lineage_name]])

        # build biom obj
        shape = [len(lineage_order), len(sample_order)],
        rows_struct = []
        cols_struct = []
        now = datetime.now()
        now.microsecond = 0
        date_str = str(now.isoformat())
        for lineage_name in lineage_order:
            rows_struct.append({'id': lineage_name, 'metadata': None})  # could add metadata full tax path if parsed from KaijuReport
        for sample_name in sample_order:
            cols_struct.append({'id': sample_name, 'metadata': None})  # sample metadata not provided to App

        biom_obj = { 'id':           output_obj_name,
                     'format':       'Biological Observation Matrix 1.0',
                     'format_url':   'http://biom-format.org/documentation/format_versions/biom-1.0.html',
                     'type':         'Taxon table',
                     'generated_by': 'KBase Kaiju App (Kaiju v1.5.0, KBase App v1.0.0)',
                     'date':         date_str,
                     'rows':         rows_struct,
                     'columns':      cols_struct,
                     'matrix_type':  'sparse',
                     'matrix_element_type': 'float',
                     'shape':        shape,
                     'data':         biom_data
                 }

        # save the biom obj to workspace
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects'] = []
        for input_reads_item in input_reads:
            provenance[0]['input_ws_objects'].append(input_reads_item['ref'])
        provenance[0]['service'] = 'kb_kaiju'
        provenance[0]['method'] = 'run_kaiju'

        if self.wsClient == None:
            try:
                self.wsClient = workspaceService(self.workspace_url, token=token)
            except:
                raise ValueError ("Unable to connect to workspace service at workspace_url: "+self.workspace_url)
        new_obj_info = wsClient.save_objects({
            'workspace':workspace_name,
            'objects':[
                { 'type': 'Communities.Biom',
                  'data': biom_obj,
                  'name': output_obj_name,
                  'meta': {},
                  'provenance': provenance
              }]
        })[0]        
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        biom_obj_ref = str(new_obj_info[WSID_I])+'/'+str(new_obj_info[OBJID_I])+'/'+str(new_obj_info[VERSION_I])

        return biom_obj_ref


    def package_folder(self, folder_path, zip_file_name, zip_file_description):
        ''' Simple utility for packaging a folder and saving to shock '''
        if folder_path == self.scratch:
            raise ValueError ("cannot package scatch itself.  folder path: "+folder_path)
        elif not folder_path.startswith(self.scratch):
            raise ValueError ("cannot package folder that is not a subfolder of scratch.  folder path: "+folder_path)
        dfu = DataFileUtil(self.callback_url)
        if not os.path.exists(folder_path):
            raise ValueError ("cannot package folder that doesn't exist: "+folder_path)
        output = dfu.file_to_shock({'file_path': folder_path,
                                    'make_handle': 0,
                                    'pack': 'zip'})
        return {'shock_id': output['shock_id'],
                'name': zip_file_name,
                'label': zip_file_description}

        
    def generate_kaijuReport_PerSamplePlots(self, options):
        pass

    def generate_kaijuReport_StackedPlots(self, options):
        tax_level = options['tax_level']
        abundance_matrix = []
        abundance_by_sample = []
        lineage_seen = dict()
        lineage_order = []
        extra_bucket_order = []
        sample_order = []
        classified_frac = []
        
        # parse summary
        for input_reads_item in options['input_reads']:
            sample_order.append(input_reads_item['name'])

            this_summary_file = os.path.join (options['in_folder'], input_reads_item['name']+'-'+tax_level+'.kaijuReport')
            (this_abundance, this_lineage_order, this_classified_frac) = self._parse_kaiju_summary_file (this_summary_file, tax_level)
            for lineage_name in this_lineage_order:
                if lineage_name not in lineage_seen:
                    lineage_seen[lineage_name] = True
                    if lineage_name.startswith('tail (<') \
                       or lineage_name.startswith('viruses') \
                       or lineage_name.startswith('unassigned at'):
                        extra_bucket_order.append(lineage_name)
                    else:
                        lineage_order.append(lineage_name)
            abundance_by_sample.append(this_abundance)
            classified_frac.append(this_classified_frac)

        for bucket_name in extra_bucket_order:  # put extra buckets at end.  necessary for sorting later.
            lineage_order.append(bucket_name)

        for lineage_i,lineage_name in enumerate(lineage_order):
            abundance_matrix.append([])
            for sample_i,sample_name in enumerate(sample_order):
                if lineage_name in abundance_by_sample[sample_i]:
                    abundance_matrix[lineage_i].append(abundance_by_sample[sample_i][lineage_name])
                else:
                    abundance_matrix[lineage_i].append(0.0)

        # make plots
        if options['plot_type'] == 'bar':
            basename_ext = '-stacked_bar_plot'
            return self._create_bar_plots (out_folder = options['stacked_plots_out_folder'], 
                                           out_file_basename = tax_level+basename_ext,
                                           vals = abundance_matrix,
                                           frac_vals = classified_frac,
                                           #title = tax_level.title()+' Level',
                                           title = tax_level.title(),
                                           frac_y_label = 'fraction classified',
                                           y_label = 'percent of classified reads',
                                           sample_labels = sample_order, 
                                           element_labels = lineage_order, 
                                           sort_by = options['sort_taxa_by'])
        elif options['plot_type'] == 'area':
            basename_ext = '-stacked_area_plot'
            return self._create_area_plots (out_folder = options['stacked_plots_out_folder'], 
                                           out_file_basename = tax_level+basename_ext,
                                           vals = abundance_matrix,
                                           frac_vals = classified_frac,
                                           #title = tax_level.title()+' Level',
                                           title = tax_level.title(),
                                           frac_y_label = 'fraction classified',
                                           y_label = 'percent of classified reads',
                                           sample_labels = sample_order, 
                                           element_labels = lineage_order, 
                                           sort_by = options['sort_taxa_by'])
        else:
            raise ValueError ("Unknown plot type "+options['plot_type'])


    def generate_kaijuReport_StackedAreaPlots(self, options):
        pass
        

    def build_html_for_kaijuReport_StackedPlots(self, input_reads, summary_folder, out_html_folder, plot_type, tax_levels, img_files):
        img_height = 750  # in pixels
        #key_scale = 25
        key_scale = img_height / 36
        img_local_path = 'img'
        out_html_img_path = os.path.join (out_html_folder, img_local_path)
        if not os.path.exists(out_html_img_path):
            os.makedirs(out_html_img_path)
        out_html_file = None
        out_html_buf = []        

        # add header
        plot_type_disp = plot_type.title()
        out_html_buf.extend (self._build_plot_html_header('KBase Kaiju Stacked '+plot_type_disp+' Abundance Plots'))

        # copy plot imgs to html folder and add img to html page
        for tax_level in tax_levels:
            src_plot_file = img_files[tax_level]
            dst_local_path = os.path.join (img_local_path, plot_type+'-'+tax_level+'.PNG')
            dst_plot_file = os.path.join (out_html_folder, dst_local_path)
            shutil.copy2 (src_plot_file, dst_plot_file)
            
            # increase height if key is long
            lineage_seen = dict()
            for input_reads_item in input_reads:
                this_summary_file = os.path.join (summary_folder, input_reads_item['name']+'-'+tax_level+'.kaijuReport')
                (this_abundance, this_lineage_order, this_classified_frac) = self._parse_kaiju_summary_file (this_summary_file, tax_level)
                for lineage_name in this_lineage_order:
                    lineage_seen[lineage_name] = True
            
            len_key = len(lineage_seen.keys())
            if key_scale * len_key > img_height:
                this_img_height = key_scale * len_key
            else:
                this_img_height = img_height

            # add img to html buf
            out_html_buf.append('<img src="'+dst_local_path+'" valign="top" height='+str(this_img_height)+'>')

        # add footer
        out_html_buf.extend (self._build_plot_html_footer())
        
        # write file
        out_local_path = plot_type+'.html'
        out_html_path = os.path.join (out_html_folder, out_local_path)
        self._write_buf_to_file(out_html_path, out_html_buf)

        out_html_file = {'type': plot_type,
                         'name': plot_type.title(),
                         'local_path': out_local_path,
                         'abs_path': out_html_path
                     }
        
        return [out_html_file]


    def build_html_for_kaijuReport_PerSamplePlots(self, out_html_folder, img_files, input_reads, tax_levels):
        img_local_path = 'img'
        out_html_img_path = os.path.join (out_html_folder, img_local_path)
        if not os.path.exists(out_html_img_path):
            os.makedirs(out_html_img_path)
        out_html_files = []

        # one page per tax_level
        for tax_level in tax_levels:
            out_html_buf = []

            # add header
            out_html_buf.extend (self._build_plot_html_header('KBase Kaiju Per-Sample Abundance Plots'))

            # copy plot imgs to html folder and add img to html page
            for input_reads_item in options['input_reads']:
                sample_name = input_reads_item['name']            
                src_plot_file = img_files[tax_level][sample_name]
                dst_local_path = os.path.join (img_local_path, 'per_sample_abundance-'+tax_level+'-'+sample_name+'.PNG')
                dst_plot_file = os.path.join (out_html_folder, dst_local_path)
                shutil.copy2 (src_plot_file, dst_plot_file)

                out_html_buf.append('<img src="'+dst_local_path+'">')

            # add footer
            out_html_buf.extend (self._build_plot_html_footer())
        
            # write file
            out_local_path = 'per_sample_abundance-'+tax_level+'.html'
            out_html_path = os.path.join (out_html_folder, out_local_path)
            self._write_buf_to_file(out_html_path, out_html_buf)

            out_html_files.append({'type': 'per_sample',
                                   'name': tax_level,
                                   'local_path': out_local_path,
                                   'abs_path': out_html_path})
                                      
        return out_html_files


    def add_top_nav(self, html_pages):
        min_downshift = 25

        for html_page in html_pages:
            html_type = html_page['type']
            name = html_page['name']
            local_path = html_page['local_path']
            abs_path = html_page['abs_path']

            # build top nav
            sp = '&nbsp;'
            sp_cnt = 2
            top_nav_buf = []
            for page_i,this_html_page in enumerate(html_pages):
                this_name = this_html_page['name']
                no_link = False
                if this_name == name:
                    no_link = True
                if this_name == 'Bar' or this_name == 'Area':
                    this_name = 'Stacked '+this_name+' Plots'
                this_local_path = this_html_page['local_path']
                if no_link:
                    disp_name = this_name.upper()
                    top_nav_item = '<b>'+disp_name+'</b>'
                else:
                    top_nav_item = '<a href="'+this_local_path+'">'+this_name+'</a>'
                if page_i == 0:
                    top_nav_item = (sp * sp_cnt) + top_nav_item
                top_nav_buf.append(top_nav_item)
            top_nav_str = ' | '.join(top_nav_buf)
            top_nav_str += '<p>'

            # add top nav to file
            new_buf = []
            with open (abs_path, 'r') as html_handle:
                for line in html_handle.readlines():
                    line_copy = line.lstrip()

                    # pad top of krona plot
                    if html_type == 'krona' and line_copy.startswith('options.style.top ='):
                        #downshift = '50px'
                        #downshift = '25px'
                        downshift_scale_per_char = 0.15
                        downshift = int(downshift_scale_per_char*len(top_nav_str))
                        if downshift < min_downshift:
                            downshift = min_downshift
                        new_buf.append("\t options.style.top = '"+str(downshift)+"px';")
                        continue

                    # capture original html
                    new_buf.append(line)

                    # add top nav str
                    if line_copy.startswith('<body'):
                        new_buf.append(top_nav_str+"\n")

            with open (abs_path, 'w') as html_handle:
                for line_buf in new_buf:
                    html_handle.write(line_buf)


    def _parse_kaiju_summary_file (self, summary_file, tax_level):
        if summary_file in self.parsed_summary:
            return (self.parsed_summary[summary_file]['abundance'],
                    self.parsed_summary[summary_file]['lineage_order'],
                    self.parsed_summary[summary_file]['classified_frac'])

        abundance = dict()
        unclassified_perc = 0.0
        unassigned_perc = None
        tail_perc = None
        virus_perc = None
        tail_cutoff = None
        lineage_order = []

        with open (summary_file, 'r') as summary_handle:
            for line in summary_handle.readlines():
                line = line.strip()
                if line.startswith('-') or line.startswith('%'):
                    continue
                (perc_str, reads_cnt_str, lineage_str) = line.split("\t")
                perc = float(perc_str.strip())
                reads_cnt = int(reads_cnt_str.strip())
                lineage = lineage_str.strip()

                if lineage == 'unclassified':
                    unclassified = perc
                elif lineage.startswith('cannot be assigned'):
                    unassigned_perc = perc
                elif lineage.startswith('belong to a'):
                    chopped_str = re.sub(r'belong to a \S+ with less than ', '', lineage)
                    tail_cutoff = re.sub(r'% of all reads', '', chopped_str)
                    tail_perc = perc
                elif lineage.startswith('Viruses'):
                    virus_perc = perc
                else:
                    lineage_order.append(lineage)
                    abundance[lineage] = perc

        if tail_cutoff != None:
            this_key = 'tail (< '+tail_cutoff+'% each taxon)'
            lineage_order.append(this_key)
            abundance[this_key] = tail_perc
        if virus_perc != None:
            this_key = 'viruses'
            lineage_order.append(this_key)
            abundance[this_key] = virus_perc
        if unassigned_perc != None:
            this_key = 'unassigned at '+tax_level+' level'
            lineage_order.append(this_key)
            abundance[this_key] = unassigned_perc

        # store to avoid repeat parse
        classified_frac = 1.0 - unclassified_perc/100.0
        self.parsed_summary[summary_file] = dict()
        self.parsed_summary[summary_file]['abundance']       = abundance
        self.parsed_summary[summary_file]['lineage_order']   = lineage_order
        self.parsed_summary[summary_file]['classified_frac'] = classified_frac

        return (abundance, lineage_order, classified_frac)


    def _create_bar_plots (self, out_folder=None, 
                           out_file_basename=None, 
                           vals=None, 
                           frac_vals=None,
                           title=None, 
                           frac_y_label=None,
                           y_label=None, 
                           sample_labels=None, 
                           element_labels=None, 
                           sort_by=None):

        # DEBUG
        #N = len(sample_labels)
        #N_force = 40
        #for sample_i in range(N_force):
        #    if sample_i >= N:
        #        sample_labels.append(sample_labels[0]+'.'+str(sample_i))
        #        frac_vals.append(frac_vals[0])
        #        for element_i,element_label in enumerate(element_labels):
        #            vals[element_i].append(vals[element_i][0])
        #for sample_i,label in enumerate(sample_labels):
        #    sample_labels[sample_i] = 'foobar.'+str(sample_i)
        # END DEBUG


        # number of samples
        N = len(sample_labels)


        # colors
        color_names = self.no_light_color_names
        len_color_names = len(color_names)
        random.seed(a=len(element_labels))
        r = random.random()
        shuffle(color_names, lambda: r)
        for label_i,label in enumerate(element_labels):
            if label_i >= len_color_names:
                color_names.append(color_names[label_i % len_color_names])
            if label.startswith('tail (<'):
                color_names[label_i] = 'lightslategray'
            elif label.startswith('viruses'):
                color_names[label_i] = 'magenta'
            elif label.startswith('unassigned at'):
                color_names[label_i] = 'darkslategray'
                

        # Sort vals
        if sort_by != None:
            print ("SORTING ELEMENTS by "+str(sort_by))
            old_index = dict()
            new_index = dict()
            for label_i,label in enumerate(element_labels):
                old_index[label] = label_i
                #print ("LABEL: "+str(label)+" OLD_INDEX: "+str(label_i))  # DEBUG

            # alphabet sort
            if sort_by == 'alpha':
                new_label_i = 0
                for label in sorted(element_labels, reverse=True):                
                    if label.startswith('tail (<') or label.startswith('viruses') or label.startswith('unassigned at'):
                        new_index[label] = old_index[label]
                    else:
                        new_index[label] = new_label_i
                        new_label_i += 1
                    #print ("LABEL: "+str(label)+" NEW_INDEX: "+str(new_index[label]))  # DEBUG

            # summed total sort
            elif sort_by == 'totals':
                totals = dict()
                for label_i,label in enumerate(element_labels):
                    totals[label] = 0
                    for sample_i,sample in enumerate(sample_labels):
                        totals[label] += vals[label_i][sample_i]
                totals_vals = []
                labels_by_totals = dict()
                for label in totals.keys():
                    if totals[label] not in totals_vals:
                        totals_vals.append(totals[label])
                        labels_by_totals[totals[label]] = []
                    labels_by_totals[totals[label]].append(label)
                new_label_i = 0
                for totals_val in sorted(totals_vals, reverse=True):
                    for label in labels_by_totals[totals_val]:
                        if label.startswith('tail (<') or label.startswith('viruses') or label.startswith('unassigned at'):
                            new_index[label] = old_index[label]
                        else:
                            new_index[label] = new_label_i
                            new_label_i += 1       
                        #print ("LABEL: "+str(label)+" NEW_INDEX: "+str(new_index[label]))  # DEBUG

            # store new order            
            new_vals = []
            new_element_labels = []
            for label_i,label in enumerate(element_labels):
                new_vals.append([])
                new_element_labels.append(None)
            for label_i,label in enumerate(element_labels):
                new_element_i = new_index[label]
                #print ("NEW_ELEMENT_I: "+str(new_element_i))  # DEBUG
                new_vals[new_element_i] = vals[label_i]
                new_element_labels[new_element_i] = label
                # DEBUG
                #print ("NEW LABEL: "+str(label)+" NEW_INDEX: "+str(new_element_i)+" OLD_INDEX: "+str(label_i))  # DEBUG
                #for sample_i,val in enumerate(new_vals[new_element_i]):
                #    print ("\t"+"SAMPLE_I: "+str(sample_i)+" NEW_VAL: "+str(new_vals[new_element_i][sample_i]))
            vals = new_vals
            element_labels = new_element_labels


        # reverse so that most important plots near top (below special 3 categories)
        element_labels = element_labels[-4::-1] + element_labels[-3:]
        vals = vals[-4::-1] + vals[-3:]

    
        # plot dimensions
        #per_unit_to_inch_scale = 0.25
        per_unit_to_inch_scale = 0.5
        bar_width_unit = 0.5
        plot_x_pad_unit = bar_width_unit / 2.0
        plot_width_unit = 2*plot_x_pad_unit + N
        downscale_above_N = 20
        extra_sample_scale = 0.5
        if N > downscale_above_N:
            plot_width_unit = 2*plot_x_pad_unit + downscale_above_N + extra_sample_scale*(N-downscale_above_N)
        plot_height_unit = 8

        
        # label dimensions
        longest_sample_label_len = 0
        longest_element_label_len = 0
        len_elements_list = len(element_labels)
        for label in sample_labels:
            if len(label) > longest_sample_label_len:
                longest_sample_label_len = len(label)
        for label in element_labels:
            if len(label) > longest_element_label_len:
                longest_element_label_len = len(label)
        #x_label_scale_unit = 0.015
        #y_label_scale_unit = 0.015
        x_label_scale_unit = 0.175
        y_label_scale_unit = 0.16
        key_label_scale = y_label_scale_unit * 50 / 30.0
        x_label_pad_unit = x_label_scale_unit * longest_element_label_len
        y_label_pad_unit = y_label_scale_unit * longest_sample_label_len
        if key_label_scale * len_elements_list > y_label_pad_unit:
            y_label_pad_unit = key_label_scale * len_elements_list
        x_label_pad_inch = per_unit_to_inch_scale * x_label_pad_unit
        y_label_pad_inch = per_unit_to_inch_scale * y_label_pad_unit


        # build canvas dimensions
        x_pad_unit = 1.0
        y_pad_unit = 0.25
        #x_pad_unit = 0.10
        #y_pad_unit = 0.10
        x_pad_inch = per_unit_to_inch_scale * x_pad_unit
        y_pad_inch = per_unit_to_inch_scale * y_pad_unit
        canvas_width_unit = 2*x_pad_unit + plot_width_unit + x_label_pad_unit
        canvas_height_unit = 2*y_pad_unit + plot_height_unit + y_label_pad_unit
        canvas_width_inch = per_unit_to_inch_scale * canvas_width_unit
        canvas_height_inch = per_unit_to_inch_scale * canvas_height_unit


        # instantiate fig
        #
        # lose axes with below grid, and so sharex property. instead match xlim, bar_width, hide ticks.
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True)  

        # gridspec_kw not in KBase docker notebook agg image (old python?)
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[1, 3]})

        # subplot2grid(shape, loc, rowspan=1, colspan=1)
        FIG_rows = 1000
        FIG_cols = 1
        top_frac = 0.22
        top_rows = int(top_frac*FIG_rows)
        bot_rows = FIG_rows-top_rows
        fig = plt.figure()
        ax_top = plt.subplot2grid((FIG_rows,FIG_cols), (0,0), rowspan=top_rows, colspan=1)
        ax_bot = plt.subplot2grid((FIG_rows,FIG_cols), (top_rows,0), rowspan=bot_rows, colspan=1)
        fig.set_size_inches(canvas_width_inch, canvas_height_inch)
        fig.tight_layout()


        #for ax in fig.axes:
        #    ax.xaxis.set_visible(False)  # remove axis labels and ticks
        #    ax.yaxis.set_visible(False)
        #    for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tick labels
        #        t.set_visible(False)
        #for ax in fig.axes:
        #    ax.spines['top'].set_visible(False) # Get rid of top axis line
        #    ax.spines['bottom'].set_visible(False) #  bottom axis line
        #    ax.spines['left'].set_visible(False) #  Get rid of bottom axis line
        #    ax.spines['right'].set_visible(False) #  Get rid of bottom axis line


        # indices
        ind = np.arange(N)    # the x locations for the groups
        label_ind = []
        for ind_i,this_ind in enumerate(ind):
            ind[ind_i] = this_ind+bar_width_unit/2
            label_ind.append(this_ind + bar_width_unit/2)
        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
        

        # plot fraction measured
        frac = ax_top.bar(ind, frac_vals, bar_width_unit, color='black', alpha=0.4, ec='none')
        ax_top.set_title(title, fontsize=11)
        ax_top.grid(b=True, axis='y')
        ax_top.set_ylabel(frac_y_label, fontsize=10)
        ax_top.tick_params(axis='y', labelsize=9, labelcolor='black')
        ax_top.set_yticks(np.arange(0.0, 1.01, .20))
        ax_top.set_ylim([0,1])
        ax_top.xaxis.set_visible(False)  # remove axis labels and ticks
        ax_top.set_xlim([-plot_x_pad_unit,N-plot_x_pad_unit])

        # plot stacked
        last_bottom = None
        p = []
        for vec_i,val_vec in enumerate(np_vals):
            if vec_i == 0:
                this_bottom = 0
                last_bottom = val_vec
            else:
                this_bottom = last_bottom
                last_bottom = this_bottom + val_vec
            p.append (ax_bot.bar (ind, val_vec, bar_width_unit, bottom=this_bottom, color=color_names[vec_i], alpha=0.4, ec='none'))

        ax_bot.set_ylabel(y_label, fontsize=10)
        ax_bot.tick_params(axis='y', direction='in', length=4, width=0.5, colors='black', labelsize=9, labelcolor='black')
        #plt.title(title)
        #plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        #ax_bot.set_xticks(label_ind, sample_labels, ha='center', rotation=90)
        ax_bot.tick_params(axis='x', direction='out', length=0, width=0, colors='black', labelsize=9, labelcolor='black')
        ax_bot.set_xticks(label_ind)
        ax_bot.set_xticklabels(sample_labels, ha='center', rotation=90)
        ax_bot.tick_params(axis='y', labelsize=9, labelcolor='black')
        ax_bot.set_yticks(np.arange(0, 101, 20))
        ax_bot.set_ylim([0,100])
        ax_bot.set_xlim([-plot_x_pad_unit,N-plot_x_pad_unit])
        

        # are positions now in units (btw. 0-1) or inches?  seems to depend on the backend Agg
        x_pad = x_pad_unit / canvas_width_unit
        y_pad = y_pad_unit / canvas_height_unit
        plot_width = plot_width_unit / canvas_width_unit
        plot_height = plot_height_unit / canvas_height_unit
        x_label_pad = x_label_pad_unit / canvas_width_unit
        y_label_pad = y_label_pad_unit / canvas_height_unit


        # Frac Plot sizing
        # don't shrink frac plot.  instead place it explictly since we built canvas for it
        #
        box = ax_top.get_position()
#        ax_top.set_position([box.x0 + x_pad_inch, 
#                             box.y0 + y_pad_inch, 
#                             box.width - x_label_pad_inch - 2*x_pad_inch, 
#                             #box.height - y_pad_inch
#                             box.height
#                         ])
        top_pos = [x_0, y_0, w, h] = [0 + x_pad, 
                                      (1.0 - top_frac)*plot_height + y_label_pad,
                                      plot_width,
                                      top_frac*plot_height - 2*y_pad
                                  ]
        ax_top.set_position(top_pos)

        # DEBUG
        #print ("AX_TOP: BOX:")
        #print ([box.x0, box.y0, box.width, box.height])
        #print ("AX_TOP: TOP_POS:")
        #print (top_pos)


        # Stacked Plot sizing
        #   don't shrink plot.  instead place it explictly since we built canvas for it
        #
        box = ax_bot.get_position()
        #ax_bot.set_position([box.x0 + x_pad_inch, 
        #                     #box.y0 + y_pad_inch + y_label_pad_inch, 
        #                     box.y0,
        #                     box.width - x_label_pad_inch - 2*x_pad_inch, 
        #                     #box.height - y_label_pad_inch - y_pad_inch
        #                     box.height
        #                 ])
        bot_pos = [x_0, y_0, w, h] = [0 + x_pad, 
                                      0 + y_label_pad + y_pad,
                                      plot_width,
                                      (1.0 - top_frac)*plot_height - 2*y_pad
                                  ]
        ax_bot.set_position(bot_pos)

        # DEBUG
        #print ("AX_BOT: BOX:")
        #print ([box.x0, box.y0, box.width, box.height])
        #print ("AX_BOT: BOT_POS:")
        #print (bot_pos)


        # add key
        key_colors = []
        for each_p in reversed(p):
            key_colors.append(each_p[0])
        ax_bot.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1,1), fontsize=9)


        # save
        img_dpi = 200
        #plt.show()
        log("SAVING STACKED BAR PLOT")
        png_file = out_file_basename+'.png'
        pdf_file = out_file_basename+'.pdf'
        output_png_file_path = os.path.join(out_folder, png_file);
        output_pdf_file_path = os.path.join(out_folder, pdf_file);
        fig.savefig(output_png_file_path, dpi=img_dpi)
        fig.savefig(output_pdf_file_path, format='pdf')
        
        return output_png_file_path


    def _create_area_plots (self, out_folder=None, 
                           out_file_basename=None, 
                           vals=None, 
                           frac_vals=None,
                           title=None, 
                           frac_y_label=None,
                           y_label=None, 
                           sample_labels=None, 
                           element_labels=None, 
                           sort_by=None):

        # number of samples
        N = len(sample_labels)


        # colors
        color_names = self.no_light_color_names
        len_color_names = len(color_names)
        random.seed(a=len(element_labels))
        r = random.random()
        shuffle(color_names, lambda: r)
        for label_i,label in enumerate(element_labels):
            if label_i >= len_color_names:
                color_names.append(color_names[label_i % len_color_names])
            if label.startswith('tail (<'):
                color_names[label_i] = 'lightslategray'
            elif label.startswith('viruses'):
                color_names[label_i] = 'magenta'
            elif label.startswith('unassigned at'):
                color_names[label_i] = 'darkslategray'
                

        # Sort vals
        if sort_by != None:
            print ("SORTING ELEMENTS by "+str(sort_by))
            old_index = dict()
            new_index = dict()
            for label_i,label in enumerate(element_labels):
                old_index[label] = label_i
                #print ("LABEL: "+str(label)+" OLD_INDEX: "+str(label_i))  # DEBUG

            # alphabet sort
            if sort_by == 'alpha':
                new_label_i = 0
                for label in sorted(element_labels, reverse=True):                
                    if label.startswith('tail (<') or label.startswith('viruses') or label.startswith('unassigned at'):
                        new_index[label] = old_index[label]
                    else:
                        new_index[label] = new_label_i
                        new_label_i += 1
                    #print ("LABEL: "+str(label)+" NEW_INDEX: "+str(new_index[label]))  # DEBUG

            # summed total sort
            elif sort_by == 'totals':
                totals = dict()
                for label_i,label in enumerate(element_labels):
                    totals[label] = 0
                    for sample_i,sample in enumerate(sample_labels):
                        totals[label] += vals[label_i][sample_i]
                totals_vals = []
                labels_by_totals = dict()
                for label in totals.keys():
                    if totals[label] not in totals_vals:
                        totals_vals.append(totals[label])
                        labels_by_totals[totals[label]] = []
                    labels_by_totals[totals[label]].append(label)
                new_label_i = 0
                for totals_val in sorted(totals_vals, reverse=True):
                    for label in labels_by_totals[totals_val]:
                        if label.startswith('tail (<') or label.startswith('viruses') or label.startswith('unassigned at'):
                            new_index[label] = old_index[label]
                        else:
                            new_index[label] = new_label_i
                            new_label_i += 1       
                        #print ("LABEL: "+str(label)+" NEW_INDEX: "+str(new_index[label]))  # DEBUG

            # store new order            
            new_vals = []
            new_element_labels = []
            for label_i,label in enumerate(element_labels):
                new_vals.append([])
                new_element_labels.append(None)
            for label_i,label in enumerate(element_labels):
                new_element_i = new_index[label]
                #print ("NEW_ELEMENT_I: "+str(new_element_i))  # DEBUG
                new_vals[new_element_i] = vals[label_i]
                new_element_labels[new_element_i] = label
                # DEBUG
                #print ("NEW LABEL: "+str(label)+" NEW_INDEX: "+str(new_element_i)+" OLD_INDEX: "+str(label_i))  # DEBUG
                #for sample_i,val in enumerate(new_vals[new_element_i]):
                #    print ("\t"+"SAMPLE_I: "+str(sample_i)+" NEW_VAL: "+str(new_vals[new_element_i][sample_i]))
            vals = new_vals
            element_labels = new_element_labels


        # reverse so that most important plots near top (below special 3 categories)
        element_labels = element_labels[-4::-1] + element_labels[-3:]
        vals = vals[-4::-1] + vals[-3:]

    
        # plot dimensions
        #per_unit_to_inch_scale = 0.25
        per_unit_to_inch_scale = 0.5
        bar_width_unit = 0.5
        plot_x_pad_unit = bar_width_unit / 2.0
        plot_width_unit = 2*plot_x_pad_unit + N
        downscale_above_N = 20
        extra_sample_scale = 0.5
        if N > downscale_above_N:
            plot_width_unit = 2*plot_x_pad_unit + downscale_above_N + extra_sample_scale*(N-downscale_above_N)
        plot_height_unit = 8

        
        # label dimensions
        longest_sample_label_len = 0
        longest_element_label_len = 0
        len_elements_list = len(element_labels)
        for label in sample_labels:
            if len(label) > longest_sample_label_len:
                longest_sample_label_len = len(label)
        for label in element_labels:
            if len(label) > longest_element_label_len:
                longest_element_label_len = len(label)
        #x_label_scale_unit = 0.015
        #y_label_scale_unit = 0.015
        x_label_scale_unit = 0.175
        y_label_scale_unit = 0.16
        key_label_scale = y_label_scale_unit * 50 / 30.0
        x_label_pad_unit = x_label_scale_unit * longest_element_label_len
        y_label_pad_unit = y_label_scale_unit * longest_sample_label_len
        if key_label_scale * len_elements_list > y_label_pad_unit:
            y_label_pad_unit = key_label_scale * len_elements_list
        x_label_pad_inch = per_unit_to_inch_scale * x_label_pad_unit
        y_label_pad_inch = per_unit_to_inch_scale * y_label_pad_unit


        # build canvas dimensions
        x_pad_unit = 1.0
        y_pad_unit = 0.25
        #x_pad_unit = 0.10
        #y_pad_unit = 0.10
        x_pad_inch = per_unit_to_inch_scale * x_pad_unit
        y_pad_inch = per_unit_to_inch_scale * y_pad_unit
        canvas_width_unit = 2*x_pad_unit + plot_width_unit + x_label_pad_unit
        canvas_height_unit = 2*y_pad_unit + plot_height_unit + y_label_pad_unit
        canvas_width_inch = per_unit_to_inch_scale * canvas_width_unit
        canvas_height_inch = per_unit_to_inch_scale * canvas_height_unit


        # instantiate fig
        #
        # lose axes with below grid, and so sharex property. instead match xlim, bar_width, hide ticks.
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True)  

        # gridspec_kw not in KBase docker notebook agg image (old python?)
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[1, 3]})

        # subplot2grid(shape, loc, rowspan=1, colspan=1)
        FIG_rows = 1000
        FIG_cols = 1
        top_frac = 0.22
        top_rows = int(top_frac*FIG_rows)
        bot_rows = FIG_rows-top_rows
        fig = plt.figure()
        ax_top = plt.subplot2grid((FIG_rows,FIG_cols), (0,0), rowspan=top_rows, colspan=1)
        ax_bot = plt.subplot2grid((FIG_rows,FIG_cols), (top_rows,0), rowspan=bot_rows, colspan=1)
        fig.set_size_inches(canvas_width_inch, canvas_height_inch)
        fig.tight_layout()


        #for ax in fig.axes:
        #    ax.xaxis.set_visible(False)  # remove axis labels and ticks
        #    ax.yaxis.set_visible(False)
        #    for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tick labels
        #        t.set_visible(False)
        #for ax in fig.axes:
        #    ax.spines['top'].set_visible(False) # Get rid of top axis line
        #    ax.spines['bottom'].set_visible(False) #  bottom axis line
        #    ax.spines['left'].set_visible(False) #  Get rid of bottom axis line
        #    ax.spines['right'].set_visible(False) #  Get rid of bottom axis line


        # indices
        ind = np.arange(N)    # the x locations for the groups
        label_ind = []
        for ind_i,this_ind in enumerate(ind):
            ind[ind_i] = this_ind+bar_width_unit/2
            label_ind.append(this_ind + bar_width_unit/2)
        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
        

        # plot fraction measured
        frac = ax_top.bar(ind, frac_vals, bar_width_unit, color='black', alpha=0.4, ec='none')
        ax_top.set_title(title, fontsize=11)
        ax_top.grid(b=True, axis='y')
        ax_top.set_ylabel(frac_y_label, fontsize=10)
        ax_top.tick_params(axis='y', labelsize=9, labelcolor='black')
        ax_top.set_yticks(np.arange(0.0, 1.01, .20))
        ax_top.set_ylim([0,1])
        ax_top.xaxis.set_visible(False)  # remove axis labels and ticks
        ax_top.set_xlim([-plot_x_pad_unit,N-plot_x_pad_unit])

        """
        ax.stackplot (ind, np_vals, colors=color_names, alpha=0.4, edgecolor='none')

        plt.ylabel(y_label)
        plt.title(title)
        plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        plt.yticks(np.arange(0, 101, 10))

        # creating the legend manually
        key_colors = []
        for color_i in reversed(np.arange(N-1)):
            key_colors.append(mpatches.Patch(color=color_names[color_i], alpha=0.4, ec='black'))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])   
        plt.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))
        """


        # plot stacked
        ax_bot.stackplot (ind, np_vals, colors=color_names, alpha=0.4, edgecolor='none')

        ax_bot.set_ylabel(y_label, fontsize=10)
        ax_bot.tick_params(axis='y', direction='in', length=4, width=0.5, colors='black', labelsize=9, labelcolor='black')
        #plt.title(title)
        #plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        #ax_bot.set_xticks(label_ind, sample_labels, ha='center', rotation=90)
        ax_bot.tick_params(axis='x', direction='out', length=0, width=0, colors='black', labelsize=9, labelcolor='black')
        ax_bot.set_xticks(label_ind)
        ax_bot.set_xticklabels(sample_labels, ha='center', rotation=90)
        ax_bot.tick_params(axis='y', labelsize=9, labelcolor='black')
        ax_bot.set_yticks(np.arange(0, 101, 20))
        ax_bot.set_ylim([0,100])
        ax_bot.set_xlim([-plot_x_pad_unit,N-plot_x_pad_unit])
        

        # are positions now in units (btw. 0-1) or inches?  seems to depend on the backend Agg
        x_pad = x_pad_unit / canvas_width_unit
        y_pad = y_pad_unit / canvas_height_unit
        plot_width = plot_width_unit / canvas_width_unit
        plot_height = plot_height_unit / canvas_height_unit
        x_label_pad = x_label_pad_unit / canvas_width_unit
        y_label_pad = y_label_pad_unit / canvas_height_unit


        # Frac Plot sizing
        # don't shrink frac plot.  instead place it explictly since we built canvas for it
        #
        box = ax_top.get_position()
#        ax_top.set_position([box.x0 + x_pad_inch, 
#                             box.y0 + y_pad_inch, 
#                             box.width - x_label_pad_inch - 2*x_pad_inch, 
#                             #box.height - y_pad_inch
#                             box.height
#                         ])
        top_pos = [x_0, y_0, w, h] = [0 + x_pad, 
                                      (1.0 - top_frac)*plot_height + y_label_pad,
                                      plot_width,
                                      top_frac*plot_height - 2*y_pad
                                  ]
        ax_top.set_position(top_pos)

        # DEBUG
        #print ("AX_TOP: BOX:")
        #print ([box.x0, box.y0, box.width, box.height])
        #print ("AX_TOP: TOP_POS:")
        #print (top_pos)


        # Stacked Plot sizing
        #   don't shrink plot.  instead place it explictly since we built canvas for it
        #
        box = ax_bot.get_position()
        #ax_bot.set_position([box.x0 + x_pad_inch, 
        #                     #box.y0 + y_pad_inch + y_label_pad_inch, 
        #                     box.y0,
        #                     box.width - x_label_pad_inch - 2*x_pad_inch, 
        #                     #box.height - y_label_pad_inch - y_pad_inch
        #                     box.height
        #                 ])
        bot_pos = [x_0, y_0, w, h] = [0 + x_pad, 
                                      0 + y_label_pad + y_pad,
                                      plot_width,
                                      (1.0 - top_frac)*plot_height - 2*y_pad
                                  ]
        ax_bot.set_position(bot_pos)

        # DEBUG
        #print ("AX_BOT: BOX:")
        #print ([box.x0, box.y0, box.width, box.height])
        #print ("AX_BOT: BOT_POS:")
        #print (bot_pos)


        """
        # add key
        key_colors = []
        for each_p in reversed(p):
            key_colors.append(each_p[0])
        ax_bot.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1,1), fontsize=9)
        """

        # create the legend manually
        w_scale = 0.8
        key_colors = []
        for color_i in reversed(np.arange(N-1)):
            key_colors.append(mpatches.Patch(color=color_names[color_i], alpha=0.4, ec='black'))
        box = ax_top.get_position()
        ax_top.set_position([box.x0, box.y0, box.width * w_scale, box.height])   
        box = ax_bot.get_position()
        ax_bot.set_position([box.x0, box.y0, box.width * w_scale, box.height])   
        ax_bot.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))


        # save
        img_dpi = 200
        #plt.show()
        log("SAVING STACKED AREA PLOT")
        png_file = out_file_basename+'.png'
        pdf_file = out_file_basename+'.pdf'
        output_png_file_path = os.path.join(out_folder, png_file);
        output_pdf_file_path = os.path.join(out_folder, pdf_file);
        fig.savefig(output_png_file_path, dpi=img_dpi)
        fig.savefig(output_pdf_file_path, format='pdf')
        
        return output_png_file_path




    def _create_area_plots_OLD (self, abundances):
        color_names = self.no_light_color_names
        
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import random
        from random import shuffle
        
        y_label = 'percent'
        title = 'Lineage Proportion'
        sample_labels = ['sample1', 'sample2', 'sample3', 'sample4', 'sample5']
        element_labels = ['OTU_1', 'OTU_2', 'OTU_3', 'OTU_4']

        N = len(sample_labels)
        random.seed(a=len(element_labels))

        r = random.random()
        shuffle(color_names, lambda: r)

        vals = [[20, 35, 20, 35, 27],
                [25, 22, 34, 20, 15],
                [45, 33, 36, 35, 48],
                [10, 10, 10, 10, 10]
            ]
        ind = np.arange(N)    # the x locations for the groups
        label_ind = ind

        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
    
        # Build image
        if N < 10:
            img_in_width = 2*N
        elif N < 20:
            img_in_width = N
        else:
            img_in_width = 20
        img_in_height = 5
        fig = plt.figure()
        fig.set_size_inches(img_in_width, img_in_height)
        ax = plt.subplot(111)

        # Let's turn off visibility of all tic labels and boxes here
        #for ax in fig.axes:
        #    ax.xaxis.set_visible(False)  # remove axis labels and tics
        #    ax.yaxis.set_visible(False)
        #    for t in ax.get_xticklabels()+ax.get_yticklabels():  # remove tics
        #        t.set_visible(False)
        #    ax.spines['top'].set_visible(False)     # Get rid of top axis line
        #    ax.spines['bottom'].set_visible(False)  # bottom axis line
        #    ax.spines['left'].set_visible(False)    # left axis line
        #    ax.spines['right'].set_visible(False)   # right axis line

        ax.stackplot (ind, np_vals, colors=color_names, alpha=0.4, edgecolor='none')

        plt.ylabel(y_label)
        plt.title(title)
        plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        plt.yticks(np.arange(0, 101, 10))

        # creating the legend manually
        key_colors = []
        for color_i in reversed(np.arange(N-1)):
            key_colors.append(mpatches.Patch(color=color_names[color_i], alpha=0.4, ec='black'))
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])   
        plt.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))

        #plt.grid()
        plt.show()


    def _build_plot_html_header(self, title):
        buf = []
        buf.append('<html>')
        buf.append('<head>')
        buf.append('<title>'+title+'</title>')

        style = '''
        <style style="text/css">
            a {
                color: #337ab7;
            }

            a:hover {
                color: #23527c;
            }

            table {
                border: 1px solid #bbb;
                border-collapse: collapse;
            }

            th, td {
                text-align: left;
                border: 1px solid #bbb;
                padding: 8px;
            }

            tr:nth-child(odd) {
                background-color: #f9f9f9;
            }

            tr:hover {
                background-color: #f5f5f5;
            }
        </style>'''

        buf.append(style)
        buf.append('</head>')
        buf.append('<body>')

        return buf


    def _build_plot_html_footer(self):
        buf = []
        buf.append('</body>')
        buf.append('</html>')
        return buf


    def _write_buf_to_file(self, filename, buf):
        with open (filename, 'w') as handle:
            for line_buf in buf:
                handle.write(line_buf+"\n")
