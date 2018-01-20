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

        log('Staged input directory: ' + input_dir)


        # 2) establish output folders
        output_dir = os.path.join(self.scratch, 'output_' + suffix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        #plots_dir = os.path.join(self.scratch, 'plot_' + suffix)
        html_dir = os.path.join(self.scratch, 'html_' + suffix)
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)

        kaiju_output_folder = os.path.join(output_dir, 'kaiju_output')
        if not os.path.exists(kaiju_output_folder):
            os.makedirs(kaiju_output_folder)

        kaijuReport_output_folder = os.path.join(output_dir, 'kaiju_report')
        if not os.path.exists(kaijuReport_output_folder):
            os.makedirs(kaijuReport_output_folder)

        kaijuReport_PerSamplePlots_output_folder = os.path.join(output_dir, 'kaiju_report_per_sample_plots')
        kaijuReport_StackedBarPlots_output_folder = os.path.join(output_dir, 'kaiju_report_stacked_bar_plots')
        kaijuReport_StackedAreaPlots_output_folder = os.path.join(output_dir, 'kaiju_report_stacked_area_plots')
        if not os.path.exists(kaijuReport_PerSamplePlots_output_folder):
            os.makedirs(kaijuReport_PerSamplePlots_output_folder)
        if not os.path.exists(kaijuReport_StackedBarPlots_output_folder):
            os.makedirs(kaijuReport_StackedBarPlots_output_folder)
        if not os.path.exists(kaijuReport_StackedAreaPlots_output_folder):
            os.makedirs(kaijuReport_StackedAreaPlots_output_folder)

        krona_output_folder = os.path.join(output_dir, 'krona_data')
        if not os.path.exists(krona_output_folder):
            os.makedirs(krona_output_folder)


        # 3) instantiate OutputBuilder
        output_folders = [ { 'name': 'kaiju_classifications',
                             'desc': 'Kaiju Classification',
                             'path': kaiju_output_folder
                         },
                           { 'name': 'kaiju_summaries',
                             'desc': 'Kaiju Summaries',
                             'path': kaijuReport_output_folder
                         },
                           { 'name': 'krona_data',
                             'desc': 'Krona Data',
                             'path': krona_output_folder
                         },
                           #{ 'name': 'per_sample_ranked_abundance_plots_PNG+PDF',
                           #  'desc': 'Per Sample Ranked Abundance Plots (PNG + PDF)',
                           #  'path': kaijuReport_PerSamplePlots_output_folder
                           #}, 
                           { 'name': 'stacked_bar_abundance_plots_PNG+PDF',
                             'desc': 'Stacked Bar Abundance Plots (PNG + PDF)',
                             'path': kaijuReport_StackedBarPlots_output_folder
                         }
                           #{ 'name': 'stacked_area_abundance_plots_PNG+PDF',
                           #  'desc': 'Stacked Area Abundance Plots (PNG + PDF)',
                           #  'path': kaijuReport_StackedAreaPlots_output_folder
                           #}
                      ]
        self.outputBuilder_client = OutputBuilder(output_folders, self.scratch, self.callback_url)


        # 4) run Kaiju in batch
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


        # 5) create Summary Reports in batch
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


        # 6) create Summary Report plots in batch
        kaijuReportPlots_options = {'input_reads':               expanded_input,
                                    'in_folder':                 kaijuReport_output_folder,
                                    'per_sample_plots_out_folder': kaijuReport_PerSamplePlots_output_folder,
                                    'stacked_bar_plots_out_folder': kaijuReport_StackedBarPlots_output_folder,
                                    'stacked_area_plots_out_folder': kaijuReport_StackedAreaPlots_output_folder,
                                    'tax_levels':                params['tax_levels']
                                    #'filter_percent':            params['filter_percent'],
                                    #'filter_unclassified':       params['filter_unclassified'],
                                    #'full_tax_path':             params['full_tax_path']
                                }
        kaijuReport_plot_files = self.run_kaijuReportPlots_batch (kaijuReportPlots_options)


        # 7) create HTML Summary Reports in batch
        #kaijuReportHTML_options = {'input_reads':               expanded_input,
        #                           'plot_files':                kaijuReport_plot_files,
        #                           'out_folder':                html_dir,
        #                           'tax_levels':                params['tax_levels']
        #}
        #self.run_kaijuReportHTML_batch (kaijuReportHTML_options)


        # 8) create Krona plots
        krona_options = {'input_reads':               expanded_input,
                         'in_folder':                 kaiju_output_folder,
                         'out_folder':                krona_output_folder,
                         'html_folder':               html_dir,
                         'db_type':                   params['db_type']
                     }
        self.run_krona_batch (krona_options)


        # 9) Package results
        output_packages = self._build_output_packages(params, self.outputBuilder_client)


        # 10) build the HTML report
        #self.outputBuilder_client.build_html_output_for_kaiju_with_krona_wf(html_dir, expanded_input)
        #if len(expanded_input) == 1:
        #    report_html_file = expanded_input[0]['name']+'.krona.html'
        #else:
        #    report_html_file = 'report.html'
        report_html_file = expanded_input[0]['name']+'.krona.html'
        html_zipped = self.outputBuilder_client.package_folder(html_dir, report_html_file, 'Kaiju abundance and Krona plots')


        # 11) save report
        report_params = {'message': '',
                         'direct_html_link_index': 0,
                         'html_links': [html_zipped],
                         'file_links': output_packages,
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


    def run_kaijuReportPlots_batch(self, options):
        input_reads = options['input_reads']
        per_sample_plot_files   = dict()
        stacked_bar_plot_files  = dict()
        stacked_area_plot_files = dict()
        for tax_level in options['tax_levels']:

#            # per sample plots
#            per_sample_plot_files[tax_level] = dict()
#            for input_reads_item in input_reads:
#                single_kaijuReportPlots_options = options
#                single_kaijuReportPlots_options['input_item'] = input_reads_item
#                single_kaijuReportPlots_options['tax_level'] = tax_level
#            
#                per_sample_plot_files[tax_level][input_reads['name']] = self.outputBuilder_client.generate_kaijuReport_PerSamplePlots(single_kaijuReportPlots_options)

            # stacked bar plots
            kaijuReportPlots_options = options
            kaijuReportPlots_options['tax_level'] = tax_level
            stacked_bar_plot_files[tax_level] = self.outputBuilder_client.generate_kaijuReport_StackedBarPlots(kaijuReportPlots_options)

#            # stacked area plots
#            kaijuReportPlots_options = options
#            kaijuReportPlots_options['tax_level'] = tax_level
#            stacked_area_plot_files[tax_level] = self.outputBuilder_client.generate_kaijuReport_StackedAreaPlots(kaijuReportPlots_options)

#        return {'per_sample_plot_files': per_sample_plot_files,
#                'stacked_bar_plot_files': stacked_bar_plot_files,
#                'stacked_area_plot_files': stacked_area_plot_files
#            }
        return {
                'stacked_bar_plot_files': stacked_bar_plot_files
            }


    def run_krona_batch(self, options, dropOutput=False):
        input_reads = options['input_reads']
        for input_reads_item in input_reads:

            # kaiju2krona
            single_kaiju2krona_run_options = options
            single_kaiju2krona_run_options['input_item'] = input_reads_item

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kaiju2krona' + '.stdout')
            
            command = self._build_kaiju2krona_command(single_kaiju2krona_run_options)
            self.run_proc (command, log_output_file)

            # kronaImport
            single_kronaImport_run_options = options
            single_kronaImport_run_options['input_item'] = input_reads_item

            log_output_file = None
            if dropOutput:  # if output is too chatty for STDOUT
                log_output_file = os.path.join(self.scratch, input_reads_item['name'] + '.kronaImport' + '.stdout')
            
            command = self._build_kronaImport_command(single_kronaImport_run_options)
            self.run_proc (command, log_output_file)


    def _validate_kaiju_options(self, options):
        # 1st order required
        func_name = 'kaiju'
        required_opts = [ 'input_item',
                          'out_folder',
                          'db_type',
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
        #options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')  # don't need for kaiju cmd
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
        required_opts = [ 'in_folder',
                          'input_item',
                          'out_folder',
                          'db_type',
                          'tax_level'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['in_folder'], options['input_item']['name']+'.kaiju')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty kaiju classification file: '+in_file)

        # db validation
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NAMES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaijuReport_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_NAMES'):
            command_list.append('-n')
            command_list.append(str(options.get('KAIJU_DB_NAMES')))
        if options.get('in_folder'):
            in_file = options['input_item']['name']+'.kaiju'
            in_path = os.path.join(options['in_folder'], in_file)
            command_list.append('-i')
            command_list.append(in_path)
        if options.get('tax_level'):
            command_list.append('-r')
            command_list.append(str(options.get('tax_level')))
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'-'+str(options.get('tax_level'))+'.kaijuReport'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        if options.get('filter_percent'):
            command_list.append('-m')
            command_list.append(str(options.get('filter_percent')))
        if int(options.get('filter_unclassified')) == 1:
            command_list.append('-u')
        if int(options.get('full_tax_path')) == 1:
            command_list.append('-p')
        

    def _build_kaijuReport_command(self, options):
        KAIJU_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_REPORT_BIN = os.path.join(KAIJU_BIN_DIR, 'kaijuReport')
        KAIJU_DB_DIR     = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')

        self._validate_kaijuReport_options(options)
        command = [KAIJU_REPORT_BIN]
        self._process_kaijuReport_options(command, options)
        return command


    def _validate_kaiju2krona_options(self, options):
        # 1st order required
        func_name = 'kaiju2krona'
        required_opts = [ 'in_folder',
                          'input_item',
                          'out_folder',
                          'db_type'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['in_folder'], options['input_item']['name']+'.kaiju')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty kaiju classification file: '+in_file)

        # db validation
        DB = 'KAIJU_DB_NODES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])
        DB = 'KAIJU_DB_NAMES'
        if not os.path.getsize(options[DB]) > 0:
            raise ValueError ('missing or empty '+DB+' file: '+options[DB])


    def _process_kaiju2krona_options(self, command_list, options):
        if options.get('KAIJU_DB_NODES'):
            command_list.append('-t')
            command_list.append(str(options.get('KAIJU_DB_NODES')))
        if options.get('KAIJU_DB_NAMES'):
            command_list.append('-n')
            command_list.append(str(options.get('KAIJU_DB_NAMES')))
        if options.get('in_folder'):
            in_file = options['input_item']['name']+'.kaiju'
            in_path = os.path.join(options['in_folder'], in_file)
            command_list.append('-i')
            command_list.append(in_path)
        if options.get('out_folder'):
            out_file = options['input_item']['name']+'.krona'
            out_path = os.path.join (str(options.get('out_folder')), out_file)
            command_list.append('-o')
            command_list.append(out_path)
        

    def _build_kaiju2krona_command(self, options):
        KAIJU_BIN_DIR     = os.path.join(os.path.sep, 'kb', 'module', 'kaiju', 'bin')
        KAIJU_2_KRONA_BIN = os.path.join(KAIJU_BIN_DIR, 'kaiju2krona')
        KAIJU_DB_DIR      = os.path.join(os.path.sep, 'data', 'kaijudb', options['db_type'])

        options['KAIJU_DB_NODES'] = os.path.join(KAIJU_DB_DIR, 'nodes.dmp')
        options['KAIJU_DB_NAMES'] = os.path.join(KAIJU_DB_DIR, 'names.dmp')

        self._validate_kaiju2krona_options(options)
        command = [KAIJU_2_KRONA_BIN]
        self._process_kaiju2krona_options(command, options)
        return command


    def _validate_kronaImport_options(self, options):
        # 1st order required
        func_name = 'kronaImport'
        required_opts = [ 'html_folder',
                          'input_item',
                          'out_folder',
                          'db_type'
                      ]
        for opt in required_opts:
            if opt not in options or options[opt] == None or options[opt] == '':
                raise ValueError ("Must define required opt: '"+opt+"' for func: '"+str(func_name)+"()'")

        # input file validation
        in_file = os.path.join(options['out_folder'], options['input_item']['name']+'.krona')
        if not os.path.getsize(in_file) > 0:
            raise ValueError ('missing or empty krona input file: '+in_file)


    def _process_kronaImport_options(self, command_list, options):
        if options.get('html_folder'):
            html_file = options['input_item']['name']+'.krona.html'
            html_path = os.path.join (str(options.get('html_folder')), html_file)
            command_list.append('-o')
            command_list.append(html_path)
        if options.get('out_folder'):
            in_file = options['input_item']['name']+'.krona'
            in_path = os.path.join(options['out_folder'], in_file)
            command_list.append(in_path)
        

    def _build_kronaImport_command(self, options):
        KRONA_BIN_DIR    = os.path.join(os.path.sep, 'usr', 'local', 'bin')
        KRONA_IMPORT_BIN = os.path.join(KRONA_BIN_DIR, 'ktImportText')
        #KRONA_BIN_DIR    = os.path.join(os.path.sep, 'kb', 'module', 'Krona', 'KronaTools', 'scripts')
        #KRONA_IMPORT_BIN = os.path.join(KRONA_BIN_DIR, 'ImportText.pl')

        self._validate_kronaImport_options(options)
        command = [KRONA_IMPORT_BIN]
        self._process_kronaImport_options(command, options)
        return command


    def _build_output_packages(self, params, outputBuilder):

        output_packages = []
        for output_folder in outputBuilder.output_folders:
            if 'skip_output_dirs' in params and output_folder['name'] in params['skip_output_dirs']:
                log('skipping output directory '+output_folder['name'])
                continue

            log('packaging output directory '+output_folder['name'])
            zipped_output_file = outputBuilder.package_folder(output_folder['path'],
                                                              output_folder['name']+'.zip',
                                                              output_folder['desc'])
            output_packages.append(zipped_output_file)

        return output_packages
