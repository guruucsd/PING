
def combine_genetic_data(export_data, gene_file):
    # Merge the two together.
    with open(gene_file, 'r') as fp:
        gene_data = np.recfromcsv(fp)

    all_idx = []
    for subj_id in export_data['SubjID']:
        x_subj_id = '"%s"' % subj_id
        if x_subj_id in gene_data['subjid']:
            all_idx.append(gene_data['subjid'].tolist().index(x_subj_id))
        else:
            all_idx.append(np.nan)

    for key in gene_data.dtype.names:
        if key == 'subjid':
            continue
        vals = [gene_data[idx] if np.logical_not(np.isnan(idx)) else np.nan
                for idx in all_idx]
        export_data[key] = vals

    return export_data
    