
from glob import glob
from src.ACESII.data_tools.data_paths import DataPaths
from src.ACESII.mission_attributes import ACESII
import spaceToolsLib as stl
import os

class DataClasses:

    def ACEII_file_executor(self, func, dict_file_path, rocket_str, just_print_file_names_bool):

        rocket_folder_path = DataPaths.ACES_data_folder
        output_file_paths = []

        # determine which rocket it is
        rocket_idx = 0 if rocket_str.lower() == 'high' else 1
        for data_name in dict_file_path.keys():
            input_path_modifier = dict_file_path[data_name][0]
            data_folder_path = f'{rocket_folder_path}/{input_path_modifier}/{data_name}/{ACESII.fliers[rocket_idx]}/'
            input_files = glob(f'{data_folder_path}/*.cdf')

            if len(input_files) == 0:
                raise Exception(f'There are no .cdf files in the {dict_file_path[data_name][0]}/{data_name}/{rocket_str} directory')
            else:
                input_names = [ifile.replace(f'{data_folder_path}/', '') for ifile in input_files]
                input_names_searchable = [ifile.replace(input_path_modifier.lower() + '_', '') for ifile in input_names]

                # Just print the file names
                if just_print_file_names_bool:
                    print(f'--- {data_name} ({rocket_str.upper()}) ---')
                    for i, file in enumerate(input_files):
                        print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
                    print('\n')

                elif dict_file_path[data_name][1][rocket_idx] == []:
                    for file_path in input_files:
                        output_file_paths.append(file_path)
                else:
                    for idx in dict_file_path[data_name][1][rocket_idx]:
                        output_file_paths.append(input_files[idx])

        if just_print_file_names_bool:
            return
        else:
            return func([stl.loadDictFromFile(file_path) for file_path in output_file_paths])


