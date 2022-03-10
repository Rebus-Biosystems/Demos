
import os
import sys 
import time
import argparse
from pathlib import Path
import pandas as pd 
import glob
from napari import Viewer, gui_qt
import random
import boto3
import botocore
from botocore import UNSIGNED
from botocore.config import Config


MIN_CRA_COL_VALUE = 2.2
MAX_CRA_COL_VALUE = 7 
DISK_SIZE = 1 

def parse_args():
    """
    Parse arguments.
    """
    # args
    parser = argparse.ArgumentParser(
        description="""
        Find dots (or dot-like objects) and write csv output file.
        """,
        epilog="""
        This pipeline will process a given experiment folder to detect RNA dots.
        """,
    )

    parser.add_argument(
        "--rootdir",
        help="The input directory containing the experiment data. ",
        metavar="dir",
        default="",
    )

    parser.add_argument(
        "--outdir",
        "-o",
        help="Specify the output directory for results files. "
        "If this is not specified then the files are written to the tmp dir.",
        metavar="dir",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    """ Main function """
    
    print("Running stitch_RNA_tables")
    stime = time.perf_counter()
    args = parse_args()

    table_files = glob.glob(os.path.join(args.rootdir, "*.csv" ))

    gene_list = [Path(table_file).stem for table_file in table_files]
    print(gene_list)

    GeneTables_folder = args.rootdir

    BUCKET_NAME = 'datafordemorebus' 
    PATH = '2021-08-13_Amsbio mouse brain_Esper HiFi_Imaging buffer_30DM genes/Tables' 

    s3 = boto3.resource('s3', config=Config(signature_version=UNSIGNED))
    bucket = s3.Bucket(BUCKET_NAME)
    # objs = bucket.objects.filter(Prefix='csv')
    for object in bucket.objects.all():
        print(object.key, object.storage_class)

    try:
        s3.Bucket(BUCKET_NAME).download_file(PATH, '*.csv')
    except botocore.exceptions.ClientError as e: 
        if e.response['Error']['Code'] == "404":
            print("The object does not exist.")
        else:
            raise



    with gui_qt():

        viewer = Viewer()
        for target in gene_list:
            print(target)
            
            target_filename = os.path.join(GeneTables_folder, target+".csv")
            if os.path.isfile(target_filename):
                results_filtered_df = pd.read_csv(target_filename)

                # TODO: Filter based on the CNR column to reduce potential FP 
                if 'CBR' in results_filtered_df.columns:
                    results_filtered_df = results_filtered_df.loc[results_filtered_df['CBR'] > MIN_CRA_COL_VALUE]
                    results_filtered_df = results_filtered_df.loc[results_filtered_df['CBR'] < MAX_CRA_COL_VALUE]
                    results_filtered_df = results_filtered_df.sort_values(by=['CBR'])

                coords_selected_df = pd.DataFrame()
                if 'x' in results_filtered_df.columns:
                    coords_selected_df['x'] = results_filtered_df['y']
                if 'y' in results_filtered_df.columns:
                    coords_selected_df['y'] = results_filtered_df['x']
                
                # Drop the dupilicates
                coords_selected_df = coords_selected_df.drop_duplicates(subset=['x','y'])

                coords_selected_arr = coords_selected_df.to_numpy()
                coords_selected = coords_selected_arr

                r = random.uniform(0,1)
                g = random.uniform(0,1)
                b = random.uniform(0.2,1)
                rgb = [r,g,b]
                layer = viewer.add_points(
                            coords_selected,
                            opacity=0.95,
                            face_color=rgb,
                            edge_color=rgb,
                            symbol="disc",
                            size=DISK_SIZE,
                            name=target,
                            )
                layer.name = target + "_CBR"
            
            else:
                print(f"{target} file does not exist in the {target_filename}" )