#!/usr/bin/env python3
"""
2D Visualization Script for Heat DNS Code
Based on ppf-main INSTOUT approach
Author: Copilot (inspired by ppf-main methodology)

This script reads both text and binary output files from the heat DNS code
and creates 2D visualizations of velocity and temperature fluctuations.

Usage:
    python visualize_heat_2d.py

The script automatically detects available data files and creates visualizations.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import struct
import os
import glob
import re

class HeatVisualization:
    """
    Class for visualizing heat DNS simulation data
    """
    
    def __init__(self, data_dir='.'):
        """
        Initialize visualization class
        
        Args:
            data_dir (str): Directory containing data files
        """
        self.data_dir = data_dir
        self.fig_dir = 'figures'
        
        # Create figures directory if it doesn't exist
        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)
            
        # Define colormaps (similar to ppf-main style)
        self.velocity_cmap = 'RdBu_r'
        self.temperature_cmap = 'coolwarm'
        
    def read_text_2d_data(self, filename, ncols, skip_header=0):
        """
        Read 2D data from text files (from TMSRSU/TMSRST output blocks 91-99)
        
        Args:
            filename (str): Name of the data file
            ncols (int): Number of columns in the file
            skip_header (int): Number of header lines to skip
            
        Returns:
            numpy.ndarray: 2D data array
        """
        try:
            data = np.loadtxt(filename, skiprows=skip_header)
            if data.ndim == 1:
                data = data.reshape(-1, ncols)
            return data
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            return None
    
    def read_binary_2d_data(self, filename):
        """
        Read binary data from BINARYOUT subroutine output
        Similar to ppf-main INSTOUT binary format
        
        Args:
            filename (str): Name of the binary data file
            
        Returns:
            dict: Dictionary containing coordinate and field data
        """
        try:
            with open(filename, 'rb') as f:
                # Read binary data (assuming single precision floats)
                data = f.read()
                
            # Unpack binary data (single precision floats)
            floats = struct.unpack(f'{len(data)//4}f', data)
            floats = np.array(floats)
            
            if 'vel' in filename:
                # Velocity data: X, Z, U, V, W, U', V', W'
                n_coords = len(floats) // 8
                ig = int(np.sqrt(n_coords))  # Assume square grid
                kg = ig
                
                # Reshape data
                x = floats[0:ig]
                z = floats[ig:ig+kg]
                u = floats[ig+kg:ig+kg+n_coords].reshape(kg, ig)
                v = floats[ig+kg+n_coords:ig+kg+2*n_coords].reshape(kg, ig)
                w = floats[ig+kg+2*n_coords:ig+kg+3*n_coords].reshape(kg, ig)
                up = floats[ig+kg+3*n_coords:ig+kg+4*n_coords].reshape(kg, ig)
                vp = floats[ig+kg+4*n_coords:ig+kg+5*n_coords].reshape(kg, ig)
                wp = floats[ig+kg+5*n_coords:ig+kg+6*n_coords].reshape(kg, ig)
                
                return {
                    'x': x, 'z': z,
                    'u': u, 'v': v, 'w': w,
                    'up': up, 'vp': vp, 'wp': wp
                }
                
            elif 'temp' in filename:
                # Temperature data: X, Z, T1, T2, T1', T2'
                n_coords = len(floats) // 6
                ig = int(np.sqrt(n_coords))
                kg = ig
                
                x = floats[0:ig]
                z = floats[ig:ig+kg]
                t1 = floats[ig+kg:ig+kg+n_coords].reshape(kg, ig)
                t2 = floats[ig+kg+n_coords:ig+kg+2*n_coords].reshape(kg, ig)
                t1p = floats[ig+kg+2*n_coords:ig+kg+3*n_coords].reshape(kg, ig)
                t2p = floats[ig+kg+3*n_coords:ig+kg+4*n_coords].reshape(kg, ig)
                
                return {
                    'x': x, 'z': z,
                    't1': t1, 't2': t2,
                    't1p': t1p, 't2p': t2p
                }
                
        except Exception as e:
            print(f"Error reading binary file {filename}: {e}")
            return None
    
    def plot_text_data_2d(self, filename_pattern, title_prefix, field_names):
        """
        Plot 2D data from text files (blocks 91-99 output)
        
        Args:
            filename_pattern (str): Pattern to match data files
            title_prefix (str): Prefix for plot titles
            field_names (list): List of field names for each column
        """
        files = glob.glob(os.path.join(self.data_dir, filename_pattern))
        files.sort()
        
        for file in files:
            # Extract time step from filename
            match = re.search(r'(\d{3})\.dat', file)
            if match:
                timestep = match.group(1)
            else:
                timestep = 'unknown'
            
            data = self.read_text_2d_data(file, len(field_names))
            if data is None:
                continue
                
            # Assume data is in format suitable for 2D grid
            # Need to determine grid dimensions - assume square grid for now
            n_points = data.shape[0]
            grid_size = int(np.sqrt(n_points))
            
            if grid_size * grid_size != n_points:
                print(f"Warning: {file} does not contain square grid data")
                continue
            
            # Create subplots for each field
            n_fields = len(field_names)
            fig, axes = plt.subplots(1, n_fields, figsize=(4*n_fields, 4))
            if n_fields == 1:
                axes = [axes]
            
            for i, (ax, field_name) in enumerate(zip(axes, field_names)):
                # Reshape data to 2D grid
                field_data = data[:, i].reshape(grid_size, grid_size)
                
                # Create 2D plot
                im = ax.imshow(field_data, 
                              extent=[0, 1, 0, 1],
                              origin='lower',
                              cmap=self.velocity_cmap if 'u' in field_name.lower() else self.temperature_cmap,
                              aspect='equal')
                
                ax.set_title(f'{field_name} (t={timestep})')
                ax.set_xlabel('x/H')
                ax.set_ylabel('z/H')
                plt.colorbar(im, ax=ax)
            
            plt.suptitle(f'{title_prefix} - Timestep {timestep}')
            plt.tight_layout()
            
            # Save figure
            basename = os.path.basename(file).replace('.dat', '')
            plt.savefig(os.path.join(self.fig_dir, f'{basename}_2d.png'), dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"Created visualization: {basename}_2d.png")
    
    def plot_binary_data_2d(self, filename_pattern, title_prefix):
        """
        Plot 2D data from binary files (BINARYOUT output)
        
        Args:
            filename_pattern (str): Pattern to match binary files
            title_prefix (str): Prefix for plot titles
        """
        files = glob.glob(os.path.join(self.data_dir, filename_pattern))
        files.sort()
        
        for file in files:
            data = self.read_binary_2d_data(file)
            if data is None:
                continue
                
            # Extract timestep from filename
            match = re.search(r'(\d+)\.dat', file)
            timestep = match.group(1) if match else 'unknown'
            
            X, Z = np.meshgrid(data['x'], data['z'])
            
            if 'vel' in file:
                # Plot velocity fluctuations
                fig, axes = plt.subplots(2, 3, figsize=(15, 10))
                
                # Mean velocities
                im1 = axes[0,0].contourf(X, Z, data['u'], levels=20, cmap=self.velocity_cmap)
                axes[0,0].set_title("Mean U")
                axes[0,0].set_aspect('equal')
                plt.colorbar(im1, ax=axes[0,0])
                
                im2 = axes[0,1].contourf(X, Z, data['v'], levels=20, cmap=self.velocity_cmap)
                axes[0,1].set_title("Mean V")
                axes[0,1].set_aspect('equal')
                plt.colorbar(im2, ax=axes[0,1])
                
                im3 = axes[0,2].contourf(X, Z, data['w'], levels=20, cmap=self.velocity_cmap)
                axes[0,2].set_title("Mean W")
                axes[0,2].set_aspect('equal')
                plt.colorbar(im3, ax=axes[0,2])
                
                # Velocity fluctuations
                im4 = axes[1,0].contourf(X, Z, data['up'], levels=20, cmap=self.velocity_cmap)
                axes[1,0].set_title("u' fluctuation")
                axes[1,0].set_aspect('equal')
                plt.colorbar(im4, ax=axes[1,0])
                
                im5 = axes[1,1].contourf(X, Z, data['vp'], levels=20, cmap=self.velocity_cmap)
                axes[1,1].set_title("v' fluctuation")
                axes[1,1].set_aspect('equal')
                plt.colorbar(im5, ax=axes[1,1])
                
                im6 = axes[1,2].contourf(X, Z, data['wp'], levels=20, cmap=self.velocity_cmap)
                axes[1,2].set_title("w' fluctuation")
                axes[1,2].set_aspect('equal')
                plt.colorbar(im6, ax=axes[1,2])
                
            elif 'temp' in file:
                # Plot temperature fluctuations
                fig, axes = plt.subplots(2, 2, figsize=(12, 10))
                
                # Mean temperatures
                im1 = axes[0,0].contourf(X, Z, data['t1'], levels=20, cmap=self.temperature_cmap)
                axes[0,0].set_title("Mean T1 (UHF)")
                axes[0,0].set_aspect('equal')
                plt.colorbar(im1, ax=axes[0,0])
                
                im2 = axes[0,1].contourf(X, Z, data['t2'], levels=20, cmap=self.temperature_cmap)
                axes[0,1].set_title("Mean T2 (CTD)")
                axes[0,1].set_aspect('equal')
                plt.colorbar(im2, ax=axes[0,1])
                
                # Temperature fluctuations
                im3 = axes[1,0].contourf(X, Z, data['t1p'], levels=20, cmap=self.temperature_cmap)
                axes[1,0].set_title("T1' fluctuation (UHF)")
                axes[1,0].set_aspect('equal')
                plt.colorbar(im3, ax=axes[1,0])
                
                im4 = axes[1,1].contourf(X, Z, data['t2p'], levels=20, cmap=self.temperature_cmap)
                axes[1,1].set_title("T2' fluctuation (CTD)")
                axes[1,1].set_aspect('equal')
                plt.colorbar(im4, ax=axes[1,1])
            
            # Set labels for all subplots
            for ax in fig.axes:
                if hasattr(ax, 'set_xlabel'):
                    ax.set_xlabel('x/H')
                    ax.set_ylabel('z/H')
            
            plt.suptitle(f'{title_prefix} - Timestep {timestep}')
            plt.tight_layout()
            
            # Save figure
            basename = os.path.basename(file).replace('.dat', '')
            plt.savefig(os.path.join(self.fig_dir, f'{basename}_contour.png'), dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"Created visualization: {basename}_contour.png")
    
    def create_all_visualizations(self):
        """
        Create all available visualizations from both text and binary data
        """
        print("Creating 2D visualizations for heat DNS simulation...")
        print(f"Data directory: {self.data_dir}")
        print(f"Output directory: {self.fig_dir}")
        
        # Text data visualizations (from blocks 91-99)
        print("\n=== Text Data Visualizations ===")
        
        # Velocity x-z plane data (block 91)
        self.plot_text_data_2d('uxz_*.dat', 'Velocity XZ Plane', ['u\' (JA)', 'u\' (JD)', 'u\' (JE)'])
        
        # Individual velocity x-z planes (blocks 92-93)
        self.plot_text_data_2d('uxz_040_*.dat', 'Velocity XZ Plane (y+≈31)', ['u\' (JD)'])
        self.plot_text_data_2d('uxz_080_*.dat', 'Velocity XZ Plane (y+≈63)', ['u\' (JE)'])
        
        # Velocity z-y plane data (block 94)
        self.plot_text_data_2d('uzy_*.dat', 'Velocity ZY Plane', ['u\'', 'v\'', 'w\''])
        
        # Velocity x-y plane data (block 95)
        self.plot_text_data_2d('uxy_*.dat', 'Velocity XY Plane', ['u\'', 'v\''])
        
        # Temperature x-z plane data (blocks 96-97)
        self.plot_text_data_2d('t1xz_*.dat', 'Temperature XZ Plane (T1&T2)', 
                              ['T1\' (JA)', 'T1\' (JD)', 'T1\' (JE)', 'T2\' (JA)', 'T2\' (JD)', 'T2\' (JE)'])
        self.plot_text_data_2d('t2xz_*.dat', 'Temperature XZ Plane (T2 only)', ['T2\' (JA)', 'T2\' (JD)', 'T2\' (JE)'])
        
        # Temperature z-y and x-y plane data (blocks 98-99)
        self.plot_text_data_2d('t12zy_*.dat', 'Temperature ZY Plane', ['T1\'', 'T2\''])
        self.plot_text_data_2d('t12xy_*.dat', 'Temperature XY Plane', ['T1\'', 'T2\''])
        
        # Binary data visualizations (from BINARYOUT)
        print("\n=== Binary Data Visualizations ===")
        self.plot_binary_data_2d('xzc_vel_*.dat', 'Velocity Field (Channel Center)')
        self.plot_binary_data_2d('xzc_temp_*.dat', 'Temperature Field (Channel Center)')
        
        print(f"\nAll visualizations saved to: {self.fig_dir}/")
        print("Visualization complete!")

def main():
    """
    Main function to run the visualization
    """
    # Initialize visualization class
    viz = HeatVisualization(data_dir='.')
    
    # Create all available visualizations
    viz.create_all_visualizations()
    
    # Show summary of created files
    print("\n=== Summary ===")
    fig_files = glob.glob(os.path.join(viz.fig_dir, '*.png'))
    print(f"Created {len(fig_files)} visualization files:")
    for f in sorted(fig_files):
        print(f"  - {os.path.basename(f)}")

if __name__ == "__main__":
    main()