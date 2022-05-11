
import os
import sys 
import time
import argparse
from pathlib import Path
import pandas as pd 
import glob
from napari import Viewer, gui_qt
import random
from cv2 import imread, cvtColor, COLOR_BGR2GRAY
import dask.array as da
import zarr 

# PARAMETERS
SHOW_FOV_GRID = False           
image_size_default = 2048
MIN_CRA_COL_VALUE = 2
MAX_CRA_COL_VALUE = 10
SHOW_NUCLEI_LOCATIONS = False 
DISK_SIZE = 1   
DISK_SIZE_NUCLEI = 4
SCALE_141 = 4
IMAGFE_SIZE_DEFAULT = 4928
DEFAULT_VISIBLE_DOTS = False 

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

    tables_folder = os.path.join(args.rootdir, "Tables") 
    stitched_image_folder = os.path.join(args.rootdir, "Stitched_images") 

    table_files = glob.glob(os.path.join(tables_folder, "*.csv" ))
    dapi_stitched = (glob.glob(os.path.join(stitched_image_folder, "nuclei_stitched_Cycle_1.tif" )))
    if len(dapi_stitched)>1:
        dapi_stitched = dapi_stitched[0]
        if not os.path.isfile(dapi_stitched):
            dapi_stitched = os.path.join(stitched_image_folder, "nuclei_stitched_coarse_cycle_1.tif")
    else:
        dapi_stitched = "None" 

    gene_list = [Path(table_file).stem for table_file in table_files]
    print(gene_list)

    GeneTables_folder = args.rootdir

    """ Visualize the RNA dots using Napari, CBR method of spot detection"""

    with gui_qt():

        viewer = Viewer()

        if os.path.isfile(dapi_stitched):
            print("Loading the DAPI image")
            nuclei_stitched = imread(dapi_stitched)
            if len(nuclei_stitched.shape)>1:
                nuclei_stitched = cvtColor(nuclei_stitched, COLOR_BGR2GRAY)

            viewer.add_image(nuclei_stitched,\
                scale=(4,4), colormap="gray_r", blending="additive", opacity=1\
                    )
        
        zarr_nuclei_stitched_filepath = os.path.join(stitched_image_folder, "nuclei_stitched_Cycle_1.zarr")
        if os.path.isdir(zarr_nuclei_stitched_filepath):
            print("Loading the DAPI Zarr image")
            nuclei_stitched_zarr = da.from_zarr(zarr_nuclei_stitched_filepath)         
            # nuclei_stitched_zarr = zarr.load(zarr_nuclei_stitched_filepath)   
            # nuclei_stitched_zarr = zarr.open(zarr_nuclei_stitched_filepath, mode='r')

            max_data_Size = 33766
            zarr_data_shap = nuclei_stitched_zarr.shape 
            y_start, y_end = 0, zarr_data_shap[1]
            x_start, x_end = 0, zarr_data_shap[0]
            if max(zarr_data_shap)>max_data_Size:
                # Divide by 4 slices 
                x_middle, y_middle = int(zarr_data_shap[0]/2), int(zarr_data_shap[1]/2)
                x_end, y_end = int(zarr_data_shap[0]), int(zarr_data_shap[1])

                for i in range(2):
                    for j in range(2):
                        viewer.add_image(
                            nuclei_stitched_zarr[i*x_middle:(i+1)*x_middle, j*y_middle:(j+1)*y_middle],
                            name="nuclei_stitched_full_res",
                            colormap="gray",
                            blending="additive",
                            translate=(i*x_middle, j*y_middle),
                            contrast_limits=(0,33000)
                        )

            else:
                viewer.add_image(
                    nuclei_stitched_zarr,
                    name="nuclei_stitched_full_res",
                    colormap="gray",
                    blending="additive",
                    contrast_limits=(0,33000)
                )
        
        for target in gene_list[0:]:
            print(target)
            
            target_filename = os.path.join(tables_folder, target + ".csv")
            if os.path.isfile(target_filename):
                results_filtered_df = pd.read_csv(target_filename)

                if 'CBR' in results_filtered_df.columns:
                    results_filtered_df = results_filtered_df.loc[results_filtered_df['CBR'] > MIN_CRA_COL_VALUE]
                    results_filtered_df = results_filtered_df.loc[results_filtered_df['CBR'] < MAX_CRA_COL_VALUE]
                    results_filtered_df = results_filtered_df.sort_values(by=['CBR'])

                coords_selected_df = pd.DataFrame()
                if 'x' in results_filtered_df.columns:
                    coords_selected_df['x'] = results_filtered_df['y']
                if 'y' in results_filtered_df.columns:
                    coords_selected_df['y'] = results_filtered_df['x']
                if 'Spot location (X)' in results_filtered_df.columns:
                    coords_selected_df['x'] = results_filtered_df['Spot location (X)']
                if 'Spot location (Y)' in results_filtered_df.columns:
                    coords_selected_df['y'] = results_filtered_df['Spot location (Y)']
                
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
                            visible = DEFAULT_VISIBLE_DOTS
                            )
                layer.name = target + "_CBR"
            
            else:
                print(f"{target} file does not exist in the {target_filename}" )

