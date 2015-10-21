"""
Add Destrieux atlas parcels (if available)
"""
import os

import numpy as np
import pandas

from .base import merge_by_key
from .default import PINGData


class DestrieuxData(PINGData):
    """ Augmented dataset with parcels based on the Destrieux (2009) atlas.

    The following files must exist in the data/Destrieux_atlas_parcels directory:
        PING_lh_area.csv, PING_rh_area.csv,
        PING_lh_thickness.csv, PING_rh_thickness.csv

    """
    IMAGING_PREFIX = PINGData.IMAGING_PREFIX + ['Destrieux_area', 'Destrieux_thickness']
    anatomical_name = dict((
        ('S_orbital-H_Shaped', 'H-shaped orbital sulcus'),
        ('G_occipital_sup', 'Superior occipital gyrus'),
        ('Medial_wall', 'medial wall (?)'),
        ('S_oc_sup_and_transversal', 'superior occipital sulcus (and transversal?)'),
        ('G_postcentral', 'postcentral gyrus'),
        ('S_calcarine', 'calcarine fissure'),
        ('G_pariet_inf-Supramar', 'inferior parietal gyrus, supramarginal'),
        ('Pole_temporal', 'temporal pole'),
        ('G_and_S_occipital_inf', 'inferior occipital sulcus/gyrus'),
        ('G_and_S_frontomargin', 'frontomarginal sulcus/gyrus'),
        ('G_precentral', 'precentral gyrus'),
        ('S_suborbital', 'suborbital sulcus'),
        ('G_and_S_transv_frontopol', 'transverse frontopolar sulcus/gyrus'),
        ('S_temporal_inf', 'inferior temporal sulcus'),
        ('S_interm_prim-Jensen', 'intermedius primus sulcus'),
        ('G_oc-temp_med-Parahip', 'medial occipito-temporal gyrus, parahippocampal'),
        ('Lat_Fis-post', 'posterior segment of the lateral fissure'),
        ('G_rectus', 'rectus gyrus'),
        ('Pole_occipital', 'occipital pole'),
        ('S_precentral-sup-part', 'precentral sulcus, superior'),
        ('S_precentral-inf-part', 'precentral sulcus, inferior'),
        ('G_front_inf-Triangul', 'IFG, pars triangularis'),
        ('S_temporal_transverse', 'transverse temporal sulcus'),
        ('G_front_sup', 'superior frontal gyrus'),
        ('Lat_Fis-ant-Vertical', 'anterior segment of lateral fissure, vertical ramus'),
        ('S_intrapariet_and_P_trans', 'intraparietal and parietal transverse sulcus'),
        ('G_and_S_cingul-Mid-Post', 'mid-posterior cingulate sulcus/gyrus'),
        ('S_parieto_occipital', 'parieto-occipital fissure'),
        ('S_temporal_sup', 'superior temporal sulcus'),
        ('S_orbital_med-olfact', 'medial orbital sulcus, olfactory'),
        ('G_and_S_subcentral', 'subcentral sulcus/gyrus'),
        ('G_temporal_inf', 'inferior temporal gyrus'),
        ('S_pericallosal', 'pericallosal sulcus'),
        ('G_subcallosal', 'subcallosal gyrus'),
        ('G_front_inf-Orbital', 'IFG, pars orbitalis'),
        ('G_front_middle', 'middle frontal gyrus'),
        ('S_postcentral', 'postcentral sulcus'),
        ('G_front_inf-Opercular', 'IFG, pars opercularis'),
        ('S_oc_middle_and_Lunatus', 'middle occipital and lunate sulcus'),
        ('S_collat_transv_ant', 'collateral transverse sulcus, anterior'),
        ('G_pariet_inf-Angular', 'inferior parietal gyrus, angular'),
        ('S_subparietal', 'subparietal sulcus'),
        ('G_temp_sup-Plan_tempo', 'planum temporale'),
        ('G_oc-temp_med-Lingual', 'medial occipital-temporal gyrus, lingual'),
        ('G_temp_sup-Plan_polar', 'planum polare'),
        ('G_orbital', 'orbital gyrus'),
        ('S_circular_insula_sup', 'circular sulcus of the insula, superior'),
        ('G_cuneus', 'cuneus'),
        ('G_parietal_sup', 'superior parietal gyrus'),
        ('G_oc-temp_lat-fusifor', 'lateral occipito-temporal gyrus'),
        ('S_cingul-Marginalis', 'cingulate sulcus, pars marginalis'),
        ('G_cingul-Post-ventral', 'post-ventral cingulate gyrus'),
        ('S_oc-temp_lat', 'occipito-temporal sulcus, lateral'),
        ('S_front_sup', 'superior frontal sulcus'),
        ('G_and_S_paracentral', 'paracentral sulcus/gyrus'),
        ('S_circular_insula_inf', 'circular sulcus of the insula, inferior'),
        ('G_occipital_middle', 'middle occipital gyrus'),
        ('S_oc-temp_med_and_Lingual', 'medial occipito-temporal sulcus, lingual'),
        ('S_collat_transv_post', 'collateral transverse sulcus, posterior'),
        ('G_temporal_middle', 'middle temporal gyrus'),
        ('G_insular_short', 'short insular gyrus'),
        ('G_precuneus', 'precuneus'),
        ('S_circular_insula_ant', 'circular sulcus of the insula, anterior'),
        ('Lat_Fis-ant-Horizont', 'anterior segment of the lateral fissure, horizontal ramus'),
        ('G_cingul-Post-dorsal', 'posterior-dorsal cingulate gyrus'),
        ('G_Ins_lg_and_S_cent_ins', 'long insular gyrus and central insular sulcus'),
        ('G_temp_sup-G_T_transv', 'superior temporal gyrus, transverse temporal gyrus'),
        ('S_occipital_ant', 'anterior occipital sulcus'),
        ('G_and_S_cingul-Mid-Ant', 'middle anterior cingulate sulcus/gyrus'),
        ('G_and_S_cingul-Ant', 'anterior cingulate sulcus/gyrus'),
        ('S_central', 'central sulcus'),
        ('S_orbital_lateral', 'lateral orbital sulcus'),
        ('S_front_middle', 'middle frontal sulcus'),
        ('G_temp_sup-Lateral', 'superior temporal gyrus, lateral'),
        ('S_front_inf', 'inferior frontal sulcus')),
        **PINGData.anatomical_name)

    anatomical_order = np.asarray(PINGData.anatomical_order.tolist() + [
        'G_and_S_transv_frontopol',
        'G_and_S_frontomargin',
        'S_orbital-H_Shaped',
        'S_orbital_lateral',
        'S_front_middle',
        'S_front_sup',
        'G_front_middle',
        'S_front_inf'
        'G_front_inf-Opercular',
        'G_front_inf-Triangul',
        'G_front_inf-Orbital',

        'S_circular_insula_ant',
        'Lat_Fis-ant-Vertical',
        'Lat_Fis-ant-Horizont',
        'S_circular_insula_sup',
        'G_insular_short',
        'G_Ins_lg_and_S_cent_ins',
        'S_circular_insula_inf',

        'G_temp_sup-Plan_polar',
        'G_temp_sup-Lateral',
        'Pole_temporal',
        'G_temporal_inf',
        'S_temporal_inf',
        'G_temporal_middle',
        'S_temporal_sup',
        'S_temporal_transverse',
        'G_temp_sup-Plan_tempo',
        'Lat_Fis-post',

        'G_pariet_inf-Supramar',
        'G_and_S_subcentral',

        'S_precentral-sup-part',
        'S_precentral-inf-part',
        'G_precentral',
        'S_central',
        'G_postcentral',
        'S_postcentral',
        'G_parietal_sup',
        'S_intrapariet_and_P_trans',
        'G_pariet_inf-Angular',
        'S_interm_prim-Jensen',
        'S_occipital_ant',

        'G_and_S_occipital_inf',
        'G_occipital_middle',
        'S_oc_sup_and_transversal',  # ?
        'S_oc_middle_and_Lunatus',
        'Pole_occipital',

        #
        'G_oc-temp_med-Parahip',
        'S_oc-temp_med_and_Lingual',
        'G_oc-temp_med-Lingual',
        'S_calcarine',
        'G_cingul-Post-ventral',
        'S_parieto_occipital',
        'G_cuneus',
        'G_occipital_sup',
        'G_precuneus',

        'S_subparietal',
        'G_cingul-Post-dorsal',
        'S_cingul-Marginalis',
        'G_and_S_cingul-Mid-Post',
        'S_pericallosal',

        'G_and_S_paracentral',
        'G_front_sup',
        # cingulate sulcus
        'G_and_S_cingul-Mid-Ant',
        'G_and_S_cingul-Ant',

        'G_orbital',
        'S_orbital_med-olfact',
        'G_subcallosal',
        'S_suborbital',
        'G_rectus',

        # ???
        'G_oc-temp_lat-fusifor',
        'G_temp_sup-G_T_transv',
        'S_oc-temp_lat',
        'S_collat_transv_ant',
        'S_collat_transv_post',
        'Medial_wall',
        ])

    def __init__(self, data=None, scrub_keys=False, scrub_values=True,
                 csv_path=None, username=None, passwd=None, force=False,
                 data_dir='data'):
        super(DestrieuxData, self).__init__(data=data, scrub_keys=scrub_keys,
                                            scrub_values=scrub_values,
                                            csv_path=csv_path, username=username,
                                            passwd=passwd, force=force, data_dir=data_dir)
        atlas_dir = os.path.join(data_dir, 'Destrieux_atlas_parcels')
        if data is not None:
            pass
        elif not os.path.exists(atlas_dir):
            raise ValueError('Destrieux data not found in %s' % atlas_dir)
        else:
            for csv_file in ['PING_lh_area.csv', 'PING_rh_area.csv',
                             'PING_lh_thickness.csv', 'PING_rh_thickness.csv']:
                csv_filepath = os.path.join(atlas_dir, csv_file)
                data = pandas.read_csv(csv_filepath)

                # Convert dots to underscores,
                # also add a prefix to each.
                measure_type = csv_file[8:-4]
                prefix = 'Destrieux_' + measure_type + '_'
                print("Converting %s data..." % csv_file)
                new_data = dict()
                for key in data:
                    val = data[key]
                    key = prefix + key.replace('_' + measure_type, '')
                    if scrub_keys:
                        key = key.replace('.', '_')
                    else:
                        key = key.replace('_rh_', '.rh.').replace('_lh_', '.lh.')
                    if scrub_values:
                        val = val.as_matrix()
                    new_data[key] = val

                # Parse and add SubjID
                id_key = [k for k in new_data if 'h.aparc.a2009s.' in k][0]
                subj_ids = [v[6:11] for v in new_data[id_key]]
                del new_data[id_key]
                new_data['SubjID'] = subj_ids
                self.merge(new_data)

    @classmethod
    def prefix2text(klass, prefix):
        prefixes = {
            'Destrieux_area.': 'Cortical surface area (mm^2)',
            'Destrieux_thickness.': 'Cortical thickness (mm)'
        }
        return prefixes.get(prefix, super(DestrieuxData, klass).prefix2text(prefix))
