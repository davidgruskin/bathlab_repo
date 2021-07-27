import numpy as np
import pandas as pd
import argparse

def main():
    # Initialize variables here:
    parser = argparse.ArgumentParser(description='Let''s clean some files!')
    parser.add_argument("-input_file", "--input_file", required=True, type=str, default='None')
    parser.add_argument("-template_file", "--template_file", required=True, type=str, default='None')
    parser.add_argument("-output_file", "--output_file", required=True, type=str, default='None')
    parser.add_argument("-outlier_threshold", "--outlier_threshold", required=True, type=float, default=.5)

    args = vars(parser.parse_args())

    # Get the filename from the first argument
    filename1 = args['input_file']
    filename2 = args['template_file']

    # Add in column names
    columns_list = ['treatment', 'sex','id', 'fam_cq_1', 'fam_cq_2', 'fam_cq_3', 'fam_sq_1', 'fam_sq_2', 'fam_sq_3',
                    'vic_cq_1', 'vic_cq_2', 'vic_cq_3', 'vic_sq_1', 'vic_sq_2', 'vic_sq_3']

    raw_data = pd.read_csv(filename1)
    template = pd.read_csv(filename2,index_col=0)

    start_idx = raw_data.index[raw_data["File Name"].str.contains("Well", na=False)].to_list()[0]
    row_dictionary = {'A': 1, 'B': 2, 'C': 3, 'D': 1, 'E': 2, 'F': 3, 'G': 1, 'H': 2, 'I': 3, 'J': 1, 'K': 2, 'L': 3,
                      'M': 1, 'N': 2, 'O': 3}
    well_counter = 0

    animal_ids = template.values.reshape(-1,).tolist()
    animal_ids = [x for x in animal_ids if str(x) != 'nan']
    animal_ids = [x for x in animal_ids if 'NTC' not in x]
    animal_ids = [x for x in animal_ids if len(x)>2]
    animal_ids = list(set(animal_ids))
    for x in range(0,len(animal_ids)):
        if len(animal_ids[x])<4:
            id_num = animal_ids[x][2]
            animal_ids[x] = animal_ids[x][0:2] + '0' + id_num
    output = pd.DataFrame(columns=columns_list,index = animal_ids)

    for row_num in range(0, raw_data.shape[0] - 19):
        well = raw_data["File Name"][row_num + 19]
        fluor = raw_data.iloc[row_num + 19, 1]
        row_id = well[0]
        replicate_num = row_dictionary[row_id]
        column_id = well[1:]
        if column_id[0] == "0":
            column_id = column_id[1]
        template_value = template.loc[row_id, str(column_id)]
        if len(template_value) < 3:
            standard_curve = 1
        elif "N" in template_value[0]:
            ntc = 1
        else:
            if len(template_value) == 3:
                id_num = template_value[2]
                template_value = template_value[0:2] + '0' + id_num
            group = template_value[0]
            sex = template_value[1]
            idd = template_value[2:]
            cq = raw_data.iloc[row_num + 19, 5]
            sq = raw_data.iloc[row_num + 19, 6]
            output.loc[template_value,"treatment"] = group
            output.loc[template_value, "sex"] = sex
            output.loc[template_value, "id"] = idd
            output.loc[template_value, fluor.lower()+"_cq"+"_"+str(replicate_num)] = cq
            output.loc[template_value, fluor.lower()+"_sq"+"_"+str(replicate_num)] = sq

    for template_value in output.index.values.tolist():
        for fluor_type in ['fam_', 'vic_']:
            for q_type in ['cq_']:
                q_value_df = output.loc[
                    template_value, [fluor_type + q_type + '1', fluor_type + q_type + '2',
                                        fluor_type + q_type + '3']].astype(float)
                if q_value_df.isna().sum() == 0:
                    max_num = float(q_value_df.max())
                    min_num = float(q_value_df.min())
                    med_num = float(q_value_df.median())
                    q_max = abs((max_num - med_num)) / abs((med_num - min_num))
                    q_min = abs((min_num - med_num)) / abs((max_num - min_num))
                    max_diff = max_num - med_num
                    min_diff = abs(min_num - med_num)
                    if max_diff > args['outlier_threshold'] and min_diff<args['outlier_threshold']:
                        if q_max > q_min:
                            idx_max = np.where(q_value_df.astype(np.float) == max_num)[0][0]
                            output.loc[template_value, [fluor_type + q_type + str(idx_max + 1)]] = '*' + output.loc[template_value, [fluor_type + q_type + str(idx_max + 1)]]
                            output.loc[template_value, [fluor_type + 'sq_' + str(idx_max + 1)]] = '*' + output.loc[template_value, [fluor_type + 'sq_' + str(idx_max + 1)]]

                    elif min_diff > args['outlier_threshold'] and max_diff< args['outlier_threshold']:
                        if q_min>q_max:
                            idx_min = np.where(q_value_df.astype(np.float) == min_num)[0][0]
                            output.loc[template_value, [fluor_type + q_type + str(idx_min + 1)]] = '*' + output.loc[template_value, [fluor_type + q_type + str(idx_min + 1)]]
                            output.loc[template_value, [fluor_type + 'sq_' + str(idx_min + 1)]] = '*' + output.loc[template_value, [fluor_type + 'sq_' + str(idx_min + 1)]]

                    elif min_diff and max_diff > args['outlier_threshold']:
                        idx_min = np.where(q_value_df.astype(np.float) == min_num)[0][0]
                        output.loc[template_value, [fluor_type + q_type + str(idx_min + 1)]] = '*' + output.loc[template_value, [fluor_type + q_type + str(idx_min + 1)]]
                        output.loc[template_value, [fluor_type + 'sq_' + str(idx_min + 1)]] = '*' + output.loc[template_value, [fluor_type + 'sq_' + str(idx_min + 1)]]

                        idx_max = np.where(q_value_df.astype(np.float) == max_num)[0][0]
                        output.loc[template_value, [fluor_type + q_type + str(idx_max + 1)]] = '*' + output.loc[template_value, [fluor_type + q_type + str(idx_max + 1)]]
                        output.loc[template_value, [fluor_type + 'sq_' + str(idx_max + 1)]] = '*' + output.loc[template_value, [fluor_type + 'sq_' + str(idx_max + 1)]]

                        idx_med = np.where(q_value_df.astype(np.float) == med_num)[0][0]
                        output.loc[template_value, [fluor_type + q_type + str(idx_med + 1)]] = '*' + output.loc[template_value, [fluor_type + q_type + str(idx_med + 1)]]
                        output.loc[template_value, [fluor_type + 'sq_' + str(idx_med + 1)]] = '*' + output.loc[template_value, [fluor_type + 'sq_' + str(idx_med + 1)]]

        output['fam_cq_average'] = output[['fam_cq_1','fam_cq_2','fam_cq_3']].apply(pd.to_numeric, args=['coerce']).mean(axis='columns',skipna=True)
        output['fam_sq_average'] = output[['fam_sq_1','fam_sq_2','fam_sq_3']].apply(pd.to_numeric, args=['coerce']).mean(axis='columns',skipna=True)
        output['vic_cq_average'] = output[['vic_cq_1','vic_cq_2','vic_cq_3']].apply(pd.to_numeric, args=['coerce']).mean(axis='columns',skipna=True)
        output['vic_sq_average'] = output[['vic_sq_1','vic_sq_2','vic_sq_3']].apply(pd.to_numeric, args=['coerce']).mean(axis='columns',skipna=True)
        output['fam_vic_sq_ratio'] = output['fam_sq_average'] /output['vic_sq_average']
        output['fam_cq_minus_vic_cq'] = output['fam_cq_average'] - output['vic_cq_average']

    output = output.sort_index()
    output.to_csv(args['output_file'], index=True)

    print("output saved to: " + args['output_file'])

if __name__ == "__main__":
    main()



