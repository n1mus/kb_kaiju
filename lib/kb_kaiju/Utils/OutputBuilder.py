import os
import shutil
import ast
import sys
import time
import re

import numpy as np
import matplotlib.pyplot as plt
import random
from random import shuffle

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

    def __init__(self, output_folders, scratch_dir, callback_url):
        self.output_folders = output_folders
        self.scratch = scratch_dir
        self.callback_url = callback_url

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

    def package_folder(self, folder_path, zip_file_name, zip_file_description):
        ''' Simple utility for packaging a folder and saving to shock '''
        if folder_path == self.scratch:
            raise ValueError ("cannot package folder that is not a subfolder of scratch")
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

    def generate_kaijuReport_StackedBarPlots(self, options):
        tax_level = options['tax_level']
        abundance_matrix = []
        abundance_by_sample = []
        lineage_seen = dict()
        lineage_order = []
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
                    lineage_order.append(lineage_name)
            abundance_by_sample.append(this_abundance)
            classified_frac.append(this_classified_frac)

        for lineage_i,lineage_name in enumerate(lineage_order):
            abundance_matrix.append([])
            for sample_i,sample_name in enumerate(sample_order):
                if lineage_name in abundance_by_sample[sample_i]:
                    abundance_matrix[lineage_i].append(abundance_by_sample[sample_i][lineage_name])
                else:
                    abundance_matrix[lineage_i].append(0.0)

        # make plots
        return self._create_bar_plots (out_folder = options['stacked_bar_plots_out_folder'], 
                                       out_file_basename = tax_level+'-stacked_bar_plot',
                                       vals = abundance_matrix,
                                       frac_vals = classified_frac,
                                       title = tax_level.title()+' Level',
                                       frac_y_label = 'fraction classified',
                                       y_label = 'percent of classified reads',
                                       sample_labels = sample_order, 
                                       element_labels = lineage_order, 
                                       sort_by = options['sort_taxa_by'])


    def generate_kaijuReport_StackedAreaPlots(self, options):
        pass
        

    def build_html_output_for_lineage_wf(self, html_dir, object_name):
        '''
        Based on the output of CheckM lineage_wf, build an HTML report
        '''

        # move plots we need into the html directory
        plot_name = 'bin_qa_plot.png'
        shutil.copy(os.path.join(self.plots_dir, plot_name), os.path.join(html_dir, plot_name))
        self._copy_ref_dist_plots(self.plots_dir, html_dir)

        # write the html report to file
        html = open(os.path.join(html_dir, 'report.html'), 'w')

        # header
        self._write_html_header(html, object_name)
        html.write('<body>\n')

        # include the single main summary figure
        html.write('<img src="' + plot_name + '" width="90%" />\n')
        html.write('<br><br><br>\n')

        # print out the info table
        self.build_summary_table(html, html_dir)

        html.write('</body>\n</html>\n')
        html.close()

        return self.package_folder(html_dir, 'report.html', 'Assembled report from CheckM')


    def build_summary_table(self, html, html_dir):

        stats_file = os.path.join(self.output_dir, 'storage', 'bin_stats_ext.tsv')
        if not os.path.isfile(stats_file):
            log('Warning! no stats file found (looking at: ' + stats_file + ')')
            return

        bin_stats = []
        with open(stats_file) as lf:
            for line in lf:
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                col = line.split('\t')
                bin_id = col[0]
                data = ast.literal_eval(col[1])
                bin_stats.append({'bid': bin_id, 'data': data})


        fields = [{'id': 'marker lineage', 'display': 'Marker Lineage'},
                  {'id': '# genomes', 'display': '# Genomes'},
                  {'id': '# markers', 'display': '# Markers'},
                  {'id': '# marker sets', 'display': '# Marker Sets'},
                  {'id': '0', 'display': '0'},
                  {'id': '1', 'display': '1'},
                  {'id': '2', 'display': '2'},
                  {'id': '3', 'display': '3'},
                  {'id': '4', 'display': '4'},
                  {'id': '5+', 'display': '5+'},
                  {'id': 'Completeness', 'display': 'Completeness', 'round': 3},
                  {'id': 'Contamination', 'display': 'Contamination', 'round': 3}]

        html.write('<table>\n')
        html.write('  <tr>\n')
        html.write('    <th><b>Bin Name</b></th>\n')
        for f in fields:
            html.write('    <th>' + f['display'] + '</th>\n')
        html.write('  </tr>\n')

        for b in bin_stats:
            html.write('  <tr>\n')
            dist_plot_file = os.path.join(html_dir, str(b['bid']) + self.DIST_PLOT_EXT)
            if os.path.isfile(dist_plot_file):
                self._write_dist_html_page(html_dir, b['bid'])
                html.write('    <td><a href="' + b['bid'] + '.html">' + b['bid'] + '</td>\n')
            else:
                html.write('    <td>' + b['bid'] + '</td>\n')
            for f in fields:
                if f['id'] in b['data']:
                    value = str(b['data'][f['id']])
                    if f.get('round'):
                        value = str(round(b['data'][f['id']], f['round']))
                    html.write('    <td>' + value + '</td>\n')
                else:
                    html.write('    <td></td>\n')
            html.write('  </tr>\n')

        html.write('</table>\n')


    # (this_abundance, this_lineage_order) = self._parse_kaiju_summary_file (this_summary_file)
    def _parse_kaiju_summary_file (self, summary_file, tax_level):
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
            this_key = 'tail < '+tail_cutoff+'% each taxon'
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

        classified_frac = 1.0 - unclassified_perc/100.0
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
        N = len(sample_labels)
        N_force = 2
        for sample_i in range(N_force):
            if sample_i >= N:
                sample_labels.append(sample_labels[0]+'.'+str(sample_i))
                frac_vals.append(frac_vals[0])
                for element_i,element_label in enumerate(element_labels):
                    vals[element_i].append(vals[element_i][0])
        # END DEBUG

        # number of samples
        N = len(sample_labels)


        # colors
        color_names = self.no_light_color_names
        random.seed(a=len(element_labels))
        r = random.random()
        shuffle(color_names, lambda: r)
        for label_i,label in enumerate(element_labels):
            if label.startswith('tail <'):
                color_names[label_i] = 'lightslategray'
            elif label.startswith('viruses'):
                color_names[label_i] = 'magenta'
            elif label.startswith('unassigned at'):
                color_names[label_i] = 'darkslategray'
                color_names = self.no_light_color_names
                

        # possibly sort vals
        if sort_by != None:
            old_index = dict()
            new_index = dict()
            for label_i,label in enumerate(element_labels):
                old_index[label] = label_i

            # alphabet sort
            if sort_by == 'alpha':
                new_label_i = 0
                for label in sorted(element_labels, reverse=True):                
                    if label.startswith('tail <') or label.startswith('viruses') or label.startswith('unassigned at'):
                        new_index[label] = old_index[label]
                    else:
                        new_index[label] = new_label_i
                        new_label_i += 1

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
                        if label.startswith('tail <') or label.startswith('viruses') or label.startswith('unassigned at'):
                            new_index[label] = old_index[label]
                        else:
                            new_index[label] = new_label_i
                            new_label_i += 1       

            # store new order            
            new_vals = []
            new_element_labels = []
            for label_i,label in enumerate(element_labels):
                new_vals.append([])
                new_element_labels.append(None)
            for label_i,label in enumerate(element_labels):
                new_vals[new_index[label]] = vals[label_i]
                new_element_labels[new_index[label]] = element_labels[label_i]
            vals = new_vals
            element_labels = new_element_labels

    
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
        for label in sample_labels:
            if len(label) > longest_sample_label_len:
                longest_sample_label_len = len(label)
        for label in element_labels:
            if len(label) > longest_element_label_len:
                longest_element_label_len = len(label)
        #x_label_scale_unit = 0.015
        #y_label_scale_unit = 0.015
        x_label_scale_unit = 0.175
        y_label_scale_unit = 0.15
        x_label_pad_unit = x_label_scale_unit * longest_element_label_len
        y_label_pad_unit = y_label_scale_unit * longest_sample_label_len
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
        top_frac = 0.15
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
                                      (1.0 - top_frac)*plot_height + y_pad + y_label_pad,
                                      plot_width,
                                      top_frac*plot_height - y_pad
                                  ]
        ax_top.set_position(top_pos)

        # DEBUG
        print ("AX_TOP: BOX:")
        print ([box.x0, box.y0, box.width, box.height])
        print ("AX_TOP: TOP_POS:")
        print (top_pos)


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
                                      (1.0 - top_frac)*plot_height - y_pad
                                  ]
        ax_bot.set_position(bot_pos)

        print ("AX_BOT: BOX:")
        print ([box.x0, box.y0, box.width, box.height])
        print ("AX_BOT: BOT_POS:")
        print (bot_pos)


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


    def _create_bar_plots_BAD_DIMENSIONS (self, out_folder=None, 
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
        N = len(sample_labels)
        N_force = 2
        for sample_i in range(N_force):
            if sample_i >= N:
                sample_labels.append(sample_labels[0]+'.'+str(sample_i))
                frac_vals.append(frac_vals[0])
                for element_i,element_label in enumerate(element_labels):
                    vals[element_i].append(vals[element_i][0])
        # END DEBUG

        # number of samples
        N = len(sample_labels)


        # colors
        color_names = self.no_light_color_names
        random.seed(a=len(element_labels))
        r = random.random()
        shuffle(color_names, lambda: r)
        for label_i,label in enumerate(element_labels):
            if label.startswith('tail <'):
                color_names[label_i] = 'lightslategray'
            elif label.startswith('viruses'):
                color_names[label_i] = 'magenta'
            elif label.startswith('unassigned at'):
                color_names[label_i] = 'darkslategray'
                color_names = self.no_light_color_names
                

        # possibly sort vals
        if sort_by != None:
            old_index = dict()
            new_index = dict()
            for label_i,label in enumerate(element_labels):
                old_index[label] = label_i

            # alphabet sort
            if sort_by == 'alpha':
                new_label_i = 0
                for label in sorted(element_labels, reverse=True):                
                    if label.startswith('tail <') or label.startswith('viruses') or label.startswith('unassigned at'):
                        new_index[label] = old_index[label]
                    else:
                        new_index[label] = new_label_i
                        new_label_i += 1

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
                        if label.startswith('tail <') or label.startswith('viruses') or label.startswith('unassigned at'):
                            new_index[label] = old_index[label]
                        else:
                            new_index[label] = new_label_i
                            new_label_i += 1       

            # store new order            
            new_vals = []
            new_element_labels = []
            for label_i,label in enumerate(element_labels):
                new_vals.append([])
                new_element_labels.append(None)
            for label_i,label in enumerate(element_labels):
                new_vals[new_index[label]] = vals[label_i]
                new_element_labels[new_index[label]] = element_labels[label_i]
            vals = new_vals
            element_labels = new_element_labels

    
        # image dimensions
        img_in_width = max_img_width = 60
        #img_in_height = max_img_height = 5
        img_in_height = max_img_height = 8
        
        if N < 5:
            img_in_width = (max_img_width/5)*N
        elif N < 10:
            img_in_width = (max_img_width/10)*N
        elif N < 20:
            img_in_width = (max_img_width/20)*N
        else:
            img_in_width = max_img_width
            img_in_height = max_img_height
        if img_in_width < 2*img_in_height:
            img_in_width = 2*img_in_height
        # scale up for later shrinkage?  No.
        #img_in_width /= 1.5*(1.0-x_shrink)
        #img_in_height /= (1.0-y_shrink)       


        # scaling based on labels
        longest_sample_label_len = 0
        longest_element_label_len = 0
        for label in sample_labels:
            if len(label) > longest_sample_label_len:
                longest_sample_label_len = len(label)
        for label in element_labels:
            if len(label) > longest_element_label_len:
                longest_element_label_len = len(label)
        max_x_shrink = 0.90
        #x_shrink_scale = 0.02*N/2.25
        #x_shrink_scale = 0.06 / N
        #x_shrink_scale = 0.01
        x_shrink_scale = 0.015
        #x_shrink_scale = 0.1
        #max_y_shrink = 0.50
        max_y_shrink = 0.20
        #y_shrink_scale = 0.0075
        y_shrink_scale = 0.015
        x_shrink = x_shrink_scale * longest_element_label_len
        y_shrink = y_shrink_scale * longest_sample_label_len
        if x_shrink > max_x_shrink:
            x_shrink = max_x_shrink
        if y_shrink > max_y_shrink:
            y_shrink = max_y_shrink


        # instantiate fig
        #
        # lose axes with below grid, and so sharex property. instead match xlim, bar_width, hide ticks.
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True)  

        # gridspec_kw not in KBase docker notebook agg image (old python?)
        #fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw = {'height_ratios':[1, 3]})

        # subplot2grid(shape, loc, rowspan=1, colspan=1)
        FIG_rows = 1000
        FIG_cols = 1
        top_frac = 0.25
        top_rows = int(top_frac*FIG_rows)
        bot_rows = FIG_rows-top_rows
        fig = plt.figure()
        ax_top = plt.subplot2grid((FIG_rows,FIG_cols), (0,0), rowspan=top_rows, colspan=1)
        ax_bot = plt.subplot2grid((FIG_rows,FIG_cols), (top_rows,0), rowspan=bot_rows, colspan=1)
        fig.set_size_inches(img_in_width, img_in_height)
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
        bar_width = 0.5      # the width of the bars: can also be len(x) sequence
        label_ind = []
        for ind_i,this_ind in enumerate(ind):
            ind[ind_i] = this_ind+bar_width/2
            label_ind.append(this_ind + bar_width/2)
        np_vals = []
        for vec_i,val_vec in enumerate(vals):
            np_vals.append(np.array(val_vec))
        

        # plot fraction measured
        frac = ax_top.bar(ind, frac_vals, bar_width, color='black', alpha=0.4, ec='none')
        ax_top.set_title(title)
        ax_top.set_ylabel(frac_y_label)
        ax_top.set_yticks(np.arange(0.0, 1.01, .20))
        ax_top.set_ylim([0,1])
        ax_top.xaxis.set_visible(False)  # remove axis labels and ticks
        ax_top.set_xlim([-0.25,N-0.25])

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
            p.append (ax_bot.bar (ind, val_vec, bar_width, bottom=this_bottom, color=color_names[vec_i], alpha=0.4, ec='none'))

        ax_bot.set_ylabel(y_label)
        #plt.title(title)
        #plt.xticks(label_ind, sample_labels, ha='right', rotation=45)
        #ax_bot.set_xticks(label_ind, sample_labels, ha='center', rotation=90)
        ax_bot.set_xticks(label_ind)
        ax_bot.set_xticklabels(sample_labels, ha='center', rotation=90)
        ax_bot.set_yticks(np.arange(0, 101, 10))
        ax_bot.set_ylim([0,100])
        ax_bot.set_xlim([-0.25,N-0.25])
        

        # Shrink frac axis
        box = ax_top.get_position()
        x_shift = 0.05
        #x_shift = 0.00
        y_shift = 0.05
        #ax_top.set_position([box.x0+x_shift, box.y0, box.width * (1.0-x_shrink), box.height-y_shift])
        ax_top.set_position([box.x0+x_shift, box.y0, box.width - x_shrink, box.height-y_shift])


        # Shrink stacked axis
        box = ax_bot.get_position()
        x_shift = 0.05
        #x_shift = 0.00
        y_shift = 0.05
        #ax_bot.set_position([box.x0+x_shift, box.y0+y_shrink, box.width * (1.0-x_shrink), box.height*(1.0-y_shrink-y_shift)])
        ax_bot.set_position([box.x0+x_shift, box.y0+y_shrink, box.width - x_shrink, box.height*(1.0-y_shrink-y_shift)])


        # add key
        key_colors = []
        for each_p in reversed(p):
            key_colors.append(each_p[0])
        ax_bot.legend(key_colors, reversed(element_labels), loc='upper left', bbox_to_anchor=(1, 1))


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


    def _create_area_plots (self, abundances):
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


    def _write_html_header(self, html, object_name):

        html.write('<html>\n')
        html.write('<head>\n')
        html.write('<title>CheckM Report for ' + object_name + '</title>')

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
        </style>\n</head>\n'''

        html.write(style)
        html.write('</head>\n')

    def _copy_file_ignore_errors(self, filename, src_folder, dest_folder):
        src = os.path.join(src_folder, filename)
        dest = os.path.join(dest_folder, filename)
        log('copying ' + src + ' to ' + dest)
        try:
            shutil.copy(src, dest)
        except:
            # TODO: add error message reporting
            log('copy failed')


    def _write_dist_html_page(self, html_dir, bin_id):

        # write the html report to file
        html = open(os.path.join(html_dir, bin_id + '.html'), 'w')

        html.write('<html>\n')
        html.write('<head>\n')
        html.write('<title>CheckM Dist Plots for Bin' + bin_id + '</title>')
        html.write('<style style="text/css">\n a { color: #337ab7; } \n a:hover { color: #23527c; }\n</style>\n')
        html.write('<body>\n')
        html.write('<br><a href="report.html">Back to summary</a><br>\n')
        html.write('<center><h2>Bin: ' + bin_id + '</h2></center>\n')
        html.write('<img src="' + bin_id + self.DIST_PLOT_EXT + '" width="90%" />\n')
        html.write('<br><br><br>\n')
        html.write('</body>\n</html>\n')
        html.close()


    def _copy_ref_dist_plots(self, plots_dir, dest_folder):
        for plotfile in os.listdir(plots_dir):
            plot_file_path = os.path.join(plots_dir, plotfile)
            if os.path.isfile(plot_file_path) and plotfile.endswith(self.DIST_PLOT_EXT):
                try:
                    shutil.copy(os.path.join(plots_dir, plotfile),
                                os.path.join(dest_folder, plotfile))
                except:
                    # TODO: add error message reporting
                    log('copy of ' + plot_file_path + ' to html directory failed')
