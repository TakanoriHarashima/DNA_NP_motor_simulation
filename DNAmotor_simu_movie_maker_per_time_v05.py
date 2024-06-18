# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:27:18 2024
@author: Takanori Harashima 
email : harashima@ims.ac.jp
"""
import os
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
from matplotlib.colors import ListedColormap
from math import log10, floor
import pandas as pd
import glob
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import subprocess
import matplotlib.patches as patches
import datetime
import time
import gc
import psutil
mplstyle.use('fast')
now = datetime.datetime.now().strftime('%Y%m%d')
start = time.time()
def round_sig(x, sig):
    try:
        return np.around(x, sig-int(floor(log10(abs(x))))-1)
    except:
        return x
def MEMORY_RELEASE():    
    gc.collect()
    mem_using = psutil.virtual_memory().percent
  #  print('Using '+ str(mem_using) + '% memory')
    if mem_using > 90:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('!!!!!!!!!!!!!       Memory error     !!!!!!!!!!!!!!!!!')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        gc.collect()
        print(memory_error_memory_error_memory_error_memory_error_memory_error)
def folder_maker(folpath):
    if os.path.isdir(folpath) == False:
        os.mkdir(folpath)
from matplotlib.collections import LineCollection
def COLOR_LINE_PLOT(x,y,c,AX): 
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='jet')
    lc.set_array(c)
    lc.set_linewidth(1.0)
    line = AX.add_collection(lc)
    return line


plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["image.interpolation"] = "none"
plt.rcParams['font.family'] = 'Arial' 
plt.rcParams['font.size'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['lines.linewidth'] = 1 
#plt.rcParams['axes.xmargin'] = 0
#plt.rcParams['axes.ymargin'] = 0

###############################################################################
RNA_size = 2000                        # Full range of RNA substrate (px)  (default:2000)
tmax_simu = 100000                     # Uplimit for simulation (sec)  (default:100000)
N_simu = 3                             # Nunber of trajectory generate (default:5)
frame_per_event = 1000                 # Span frame to check the progress of simulation (default:1000)
####################################################################################
####################################################################################
# Geometric parameters (Fixed parameters)
particle_size = 100                    # Diameter of AuNP [nm]
r_particle = particle_size / 2         # Radius of AuNP [nm]
RNA_density =0.022                     # Surface density of RNAs on substrate[molecules/nm2]
DNA_density =0.10                      # Surface density of DNAs on AuNP [molecules/nm2]
L_DNA = 7.3                            # 45 base ssDNA, end-to-end length [nm]
R_mobile_nm = 25.2                     # Mobile region of DNA/RNA duplex [nm]                
                                       # Assuming only the rotation of AuNP
                                       # Mobile region of the particle centroid when DNA/RNA duplex can be elongated to its contour length.
# simulation parameter calculation
# Accessible radius [nm]
path_radius = ((L_DNA + r_particle)**2 - r_particle**2 )**0.5 
# Averaged number of accessible RNA molecules  
path_RNA_molecules = RNA_density * (np.pi*path_radius **2)
# Accessible AuNP surface area
path_DNA_radius_nm2 = (4*np.pi*r_particle**2) * (np.degrees(np.arcsin(path_radius/(L_DNA + r_particle)))*2/360)
# Averaged number of accessible DNA molecules   = Accessible AuNP surface area x DNA surface density
path_DNA_molecules = DNA_density * path_DNA_radius_nm2
# Pixel size = 1/DNA surface density
pixel_to_nm = np.sqrt(1/DNA_density) # nm
# mobile region [pixel]
R_mobile_pixel = R_mobile_nm / pixel_to_nm
# Accessible radius [pixel]
r_path = path_radius / pixel_to_nm  
DNA_size = round(r_path) 
###############################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################# variable parameters ############################# 
############################# variable parameters ############################# 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################
workfol = r'Directory\to\perform\DNAmotor\simulation\001_khyb=0.30_kcatE=4.0_konE=1.0x106'
path_list = [workfol + os.sep + r"DATE_{'N_simu'= 3, 'tmax'= 100000, 'RNaseH'= 36, 'frame_per_event'= 1000}"]
fps = 20                         # Frame per second, adjust to experimental conditions
dt_timecourse = 500              # Range to show time-course [sec]
itrace= 0                        # Trajectory ID to make the movie 
MAKE_FIGURES = True # If True, the program generate every snapshots during simulation as png data.  
MAKE_MOVIE = True   # If True, the program concatenate the images by ffmpeg.  
###############################################################################
###############################################################################
###############################################################################
###############################################################################
dt = 1/fps * 10
colortable_RNA = ListedColormap(['0.7','g','r','w']) # ['Unbound', 'Bound', 'Hydrolyzed', 'Empty']

if MAKE_FIGURES:
    for path in path_list:
        print(path)
        datas = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*.csv'))
        datas =  [f for f in datas if f[-8:-4].isnumeric()]
        datas_RNA = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*_RNAstate.gz'))
        datas = [datas[itrace]]
        datas_RNA = [datas_RNA[itrace]]
        
        for data, data_RNA in zip(datas, datas_RNA):
            print( data, data_RNA)
            
            if os.path.exists(data.replace('.csv','.mp4')):
                continue
            
            pi = data[-8:-4]
            trace_fol = path + os.sep + 'progress' +  os.sep + pi
            folder_maker(trace_fol)
    
            df = pd.read_csv(data)
           # RNA_state = np.array(pd.read_csv(data_RNA, index_col = 0))
            RNA_state = np.array(pd.read_pickle(data_RNA))

            
            RNA_state = np.where(RNA_state != 3, 0, 3)
            
            df_param = pd.read_csv(path + os.sep + 'parameters.csv')
            
            RNaseH = int(df_param['RNaseH'])
            cx, cy = len(RNA_state)//2, len(RNA_state)//2
    
            display_size = DNA_size*3
            min_disp, max_disp =  cx-display_size, cx+display_size
            
            fi = 0
            tmin,tmax = 0, dt_timecourse
            tmin_event = list(df[df['t'] <= tmin]['t'])[-1]
            tmax_event = list(df[df['t'] >= min(max(df['t']),tmax)]['t'])[0]
            df_show = df[(df['t'] >= tmin_event) & (df['t'] <= tmax_event)]
            
            for ti in np.arange(0,max(df['t']), dt):
                
                if ti < dt_timecourse/2:
                    tmin_event = 0
                    tmax_event = tmax           
                else:
                    tmin_event = ti - dt_timecourse
                    tmax_event = ti + dt_timecourse  
                df_show = df[(df['t'] >= tmin_event) & (df['t'] <= tmax_event)]
                
                df_i = df[df['t'] <= ti]
                ei = len(df_i)
                cx_list =  np.array(df_i['X (pixel)'])
                cy_list =  np.array(df_i['Y (pixel)'])
                cx = cx_list[-1]
                cy = cy_list[-1]


                df_event_i = df_show[(df_show['t'] >= ti-dt) & (df_show['t'] <= ti)]
                for index, row in df_event_i.iterrows():
                    event_position = tuple([int(row['ex']), int(row['ey'])])
                    if row['event'] == 1:                
                        RNA_state[event_position] = 1
                    if row['event'] == 3:                
                        RNA_state[event_position] = 2
    
    
                xmin_disp, xmax_disp= min([cx-display_size,min(cx_list) ]), max([cx+display_size,max(cx_list) ])
                ymin_disp, ymax_disp= min([cy-display_size,min(cy_list) ]), max([cy+display_size,max(cy_list) ])
                min_disp, max_disp = min([xmin_disp,ymin_disp, min_disp]), max([xmax_disp,ymax_disp, max_disp])
                min_disp_nm, max_disp_nm = (min_disp-RNA_size//2)*pixel_to_nm, (max_disp-RNA_size//2)*pixel_to_nm


                fig, ax = plt.subplots(figsize=(18, 11))
                gs = GridSpec(3, 3, left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.35, hspace=0.05)
                AX1 = plt.subplot(gs[:2, 0])
                AX2 = plt.subplot(gs[:2, 1])
                AX3 = plt.subplot(gs[:2, 2])
                AX4 = plt.subplot(gs[2:, :])
                
                """                
                for AXi in [AX2, AX3]:                      
                    AXi.text(min_disp_nm+(max_disp_nm-min_disp_nm)*0.55 ,max_disp_nm-(max_disp_nm-min_disp_nm)*0.85, 'Event = '+str(ei)+'\n'+'t = '+ str(round_sig(ti,3)) + ' s' ,color = 'k')
                """
                
                cx_nm = (cx-RNA_size/2)*pixel_to_nm 
                cy_nm = (cy-RNA_size/2)*pixel_to_nm 
                cx_nm_list = (np.array(cx_list)-RNA_size/2) * pixel_to_nm 
                cy_nm_list = (np.array(cy_list)-RNA_size/2) * pixel_to_nm 

                #######################################################################################################################
                #######################################################################################################################
           #     AX1.set_title("p" + pi + ', [RNaseH]='+ str(RNaseH) + 'nM'+'\n')
                zoom_size = 20
                zoom_size_nm = zoom_size * pixel_to_nm
                #AX1.text(cx_nm+zoom_size_nm*0.20 ,cy_nm-zoom_size_nm*0.7, 'Event = '+str(ei)+'\n'+'t = '+ str(round_sig(ti,3)) + ' s' ,color = 'k', ha='left')
                
                scalebar = 20 # nm
                AX1.plot( [cx_nm+zoom_size_nm*0.5,cx_nm+zoom_size_nm*0.5 + scalebar],[cy_nm+zoom_size_nm*0.7,cy_nm+zoom_size_nm*0.7], lw = 8,c='k')
                AX1.text(cx_nm+zoom_size_nm*0.5+ scalebar/2 ,cy_nm+zoom_size_nm*0.83, str(scalebar)+' nm', va='center', ha='center',weight='bold',fontsize= 20,zorder=10)
                m1 = AX1.imshow(RNA_state[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, clim=(0,4),cmap=colortable_RNA\
                           ,extent = ( cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm, cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm), origin='lower')    
                AX1.set(xlim=(cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm),ylim=(cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm) )

                # set color bar
                AX1_inset = inset_axes(AX1,
                                    width="6%",  # width = 50% of parent_bbox width
                                    height="30%",  # height : 5%
                                    bbox_to_anchor=(0.6, 1.01, 1, 1),
                                    bbox_transform=AX1.transAxes,
                                    loc='lower left')
                cbar = fig.colorbar(m1, cax=AX1_inset, ticks=[0.5,1.5,2.5,3.5])
                cbar.ax.set_yticklabels(['Unbound', 'Bound', 'Hydrolyzed', 'Empty'])
                #AX1_inset.xaxis.set_ticks_position("top")
                #######################################################################################################################
                #######################################################################################################################

                m2 = AX2.imshow(RNA_state[min_disp:max_disp,min_disp:max_disp].T, clim=(0,4),cmap=colortable_RNA\
                           ,extent = ( min_disp_nm, max_disp_nm, min_disp_nm, max_disp_nm), origin='lower')    


                AX3.plot(cx_nm_list,cy_nm_list, lw = 1,color = '0.2',marker= 'o',ms = 2)
                AX3.plot(cx_nm,cy_nm, lw = 1,color = 'y',marker= '+',ms = 10)
                AX3.set_xlim( min_disp_nm, max_disp_nm)
                AX3.set_ylim( min_disp_nm, max_disp_nm)
                AX3.set_aspect('equal')
                
                for AXi in [AX1,AX2,AX3]:
                    AXi.axhline(0, lw = 1,color = 'gray',zorder=1)
                    AXi.axvline(0, lw = 1,color = 'gray',zorder=1)
                    c = patches.Circle(xy=( cx_nm, cy_nm), radius=r_path*pixel_to_nm, ec='y',fill=False,lw=2)
                    AXi.add_patch(c)
                    AXi.invert_yaxis()
                    AXi.set(xlabel= 'X (nm)',ylabel= 'Y (nm)')
                    
                #AX2.set_title('Event = '+str(ei)+'\n'+'t = '+ str(round(ti,2)) + ' s' ,color = 'k',loc='left',fontsize=20,va ='bottom')
                AX1.set_title('t = '+ str(round(ti,2)) + ' s' ,color = 'k',loc='left',fontsize=24,va ='bottom')
                
                AX4.axvline(x=ti, ls = 'dashed',c='0.5', lw = 1)
                AX4.plot([ti], [74], marker= "v",c='0.5',ms=12, clip_on=False)
                
                AX4.plot([ti], [df_i['N1'].iloc[-1]], marker= "o",c='0.2',ms=10)
                AX4.plot([ti], [df_i['N2'].iloc[-1] + df_i['N3'].iloc[-1]], marker= "o",c='g',ms=10)
                AX4.plot([ti], [df_i['N4'].iloc[-1]], marker= "o",c='r',ms=10)
                
                AX4.step(df_show['t'], df_show['N1'], where='post',color = '0.2',label ='Unbound',lw=2)
                AX4.step(df_show['t'], df_show['N2'] + df_show['N3'], where='post', color = 'g',label ='Bound',lw=2)
                AX4.step(df_show['t'], df_show['N4'],where='post',color = 'r',label ='Hydrolyzed',lw=2)

                if ti < dt_timecourse/2:
                    AX4.set(ylabel='Number of sites',xlabel='Time (s)',ylim=(-2,70),yticks=[0,10,20,30,40,50,60,70],xlim=(tmin,tmax) )            
                else:
                    AX4.set(ylabel='Number of sites',xlabel='Time (s)',ylim=(-2,70),yticks=[0,10,20,30,40,50,60,70],xlim=(ti-dt_timecourse/2,ti+dt_timecourse/2) )            
                
                AX4.grid(which = 'both',axis='y')
                AX4.legend(loc ='upper right', ncol= 3)                    
                #plt.subplots_adjust(hspace=0.05, wspace=0.05) 
                konDNA = round(df_param['khyb'].iloc[0],1)
                print(konDNA)
                konE   = round(df_param['konE'].iloc[0] * 10**-6,1)
                kcatE  = round(df_param['kcat'].iloc[0],1)
                title_condition = '$k_{\mathrm{on}}^\mathrm{E}$'+f' = {konE}'+r' $\times$ 10$^6$ (M$^{-1}$ s$^{-1}$)'+ '\n' + \
                             '$k_{\mathrm{cat}}^\mathrm{E}$' + f' = {kcatE}' + r' (s$^{-1}$)'+ '\n' + \
                             '$k_{\mathrm{on}}^\mathrm{DNA/RNA}$' + f' = {konDNA}' + r' (s$^{-1}$)'+ '\n' + \
                             '[RNase H] = ' +str(int(RNaseH)) + ' nM'
               # fig.suptitle(title_condition,fontsize= 24)
                AX2.set_title(title_condition,fontsize= 24,loc = 'left',va = 'bottom')
                
                plt.savefig(trace_fol + os.sep + pi +'_'+ ('00000'+str(fi))[-6:] +'.png' ,dpi = 100)

                if fi%50 == 0:
                    print('############################################################')
                    #print(data)
                    print ('i='+str(ei), ' t='+str(ti))
                    plt.show()
                    MEMORY_RELEASE()
                
                plt.clf()
                plt.close()
                ###################
                ###################
                fi += 1
                ###################
                ###################

    
if MAKE_MOVIE:
    for path in path_list:
        datas = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*.csv'))
        datas =  [f for f in datas if f[-8:-4].isnumeric()]
    
    
        for data in datas:
            pi = data[-8:-4]
            
            command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -q 8 -crf 0 '+pi+ '.avi' 
            command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -pix_fmt yuv420p '+pi+ '.mp4' 
            print(command)
            p = subprocess.run(command, stdout=subprocess.PIPE, shell=True, cwd=path + os.sep + 'progress',universal_newlines=True)
        
        
    
