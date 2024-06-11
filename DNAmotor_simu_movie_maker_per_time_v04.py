# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:27:18 2024

@author: Takanori Harashima 
email : harashima@ims.ac.jp
"""

import itertools
import tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
import cv2
from PIL import Image
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
import matplotlib.style as mplstyle
mplstyle.use('fast')
import math
from math import log10, floor
import pandas as pd
import glob
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D   

from tqdm import tqdm
import subprocess
import sys
import random
import matplotlib.animation as animation
import matplotlib.patches as patches
import cv2
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import LinearSegmentedColormap
import datetime
now = datetime.datetime.now().strftime('%Y%m%d')
import time
import gc
import psutil
start = time.time()
def round_sig(x, sig):
    try:
        return np.around(x, sig-int(floor(log10(abs(x))))-1)
    except:
        return x
def MEMORY_RELEASE():    
    #メモリ解放処理
    gc.collect()
    #メモリ使用量を確認　→　90%以上だったらプログラムを強制終了させる
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
  #  x    = np.linspace(0,1, 100)
 #   y    = np.linspace(0,1, 100)
  #  cols = np.linspace(0,1,len(x))    
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
plt.rcParams['font.family'] = 'Arial' #全体のフォントを設定
plt.rcParams['font.size'] = 16 #フォントサイズを設定
plt.rcParams['axes.titlesize'] = 16 #フォントサイズを設定
plt.rcParams['lines.linewidth'] = 1 
#plt.rcParams['axes.xmargin'] = 0
#plt.rcParams['axes.ymargin'] = 0


###############################################################################
RNA_size = 3000 # シミュレーションに使うRNA表面の範囲
tmax_simu = 100000 # sec シミュレーションする上限時間 
N_simu = 10 # simulationするサンプル数
frame_per_event = 1000 # v11 gifアニメの1フレームあたりのイベント数
####################################################################################
####################################################################################
#金ナノparticle直径
particle_size = 100 # nm
r_particle = particle_size / 2 # nm  
# RNA被覆率 Supplementary Figure 4.　https://onlinelibrary.wiley.com/doi/10.1002/anie.201916281
RNA_density =0.022 # molecules/nm2
# DNA被覆率 Figure S6. https://pubs.acs.org/doi/10.1021/acsnano.0c10658
DNA_density =0.10 # molecules/nm2

# DNA legのアクセス可能領域
L_DNA = 7.3 # nm # 45 base, end-to-end length
# DNARNAhybridの可動域
R_mobile_nm = 25.2 # nm # 粒子の転がり運動のみを仮定、真下のDNA/RNAのcontour lengthまで可動とした際の粒子重心の可動半径

####################################################################################
####################################################################################
# simulation parameter 計算
# 反応半径（nm）
path_radius = ((L_DNA + r_particle)**2 - r_particle**2 )**0.5 
# 反応可能なRNAの分子数
path_RNA_molecules = RNA_density * (np.pi*path_radius **2)
# 反応可能な「粒子上の」面積 = 4πr^2 * (カバーする角度範囲)
path_DNA_radius_nm2 = (4*np.pi*r_particle**2) * (np.degrees(np.arcsin(path_radius/(L_DNA + r_particle)))*2/360)
# 反応可能なDNA分子数 = 反応可能な「粒子上の」面積 x DNA密度
path_DNA_molecules = DNA_density * path_DNA_radius_nm2
# pixel size = 1/DNA密度, 円の真下部分の密度を用いて反応エリア内全体のDNAの密度として近似する
pixel_to_nm = np.sqrt(1/DNA_density) # nm
# RNA/DNA二重鎖のmobile region
R_mobile_pixel = R_mobile_nm / pixel_to_nm
# 反応半径
r_path = path_radius / pixel_to_nm  
DNA_size = round(r_path) 

print ("path_DNA_molecules",path_DNA_molecules)
print ("path_RNA_molecules",path_RNA_molecules)
print ('R_mobile_pixel', R_mobile_pixel)
print ('R_mobile nm', R_mobile_nm)

print('path_radius (nm)',path_radius)
print ('r_path nm', r_path *pixel_to_nm)
print ('r_path px', r_path )
print ("DNA_size",DNA_size)
####################################################################################
###############################################################################
# simulation v5.01以降
# RNAの分布, eventの起きた座標,eventの種類, 重心の座標, さえわかれば動画は作製できる
# シミュレーションの結果のcsvファイルからmovieを作成する
# v03 eventごとではなく単位時間でmovieを作成する
colortable_RNA = ListedColormap(['0.7','g','r','w']) # ['Unbound', 'Bound', 'Hydrolyzed', 'Empty']
# *RNA_state.csv : 最終フレームののデータ → 0以外を1に初期化する
globalfol = r'G:\マイドライブ\DNAmotor_simulation\DNAmotor_simu_v5.01\240129_movie' # この直下のファイルを選ぶ
workfol = globalfol + os.sep + '001_khyb=0.30_kcatE=4.0_konE=1.0x106'
path_list = sorted(glob.glob(workfol + '/*'),key=len)
print(path_list)

path_list = [r"Example_khyb=0.30_kcatE=4.0_konE=1.0x106\20240123_{'N_simu'= 20, 'tmax'= 100000, 'RNaseH'= 36, 'frame_per_event'= 1000}"]
fps = 20
dt = 1/20 * 10
dt_timecourse = 500 # s
itrace= 3
###############################################################################
###############################################################################
MAKE_FIGURES = True
MAKE_MOVIE = False
###############################################################################
###############################################################################

if MAKE_FIGURES:
    for path in path_list:
        print(path)
        datas = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*.csv'))
        datas =  [f for f in datas if f[-8:-4].isnumeric()]
        datas_RNA = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*_RNAstate.csv'))
        
        datas = [datas[itrace]]
        datas_RNA = [datas_RNA[itrace]]
        
        for data, data_RNA in zip(datas, datas_RNA):
            print( data, data_RNA)
            
            """
            # 動画がすでにあればSkip
            if os.path.exists(data.replace('.csv','.mp4')):
                continue
            """
            
            pi = data[-8:-4]
            trace_fol = path + os.sep + 'progress' +  os.sep + pi
            folder_maker(trace_fol)
    
            df = pd.read_csv(data)
            RNA_state = np.array(pd.read_csv(data_RNA, index_col = 0))
            
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
                
                """
                if tmin != ti//dt_timecourse*dt_timecourse:
                    tmin,tmax = ti//dt_timecourse*dt_timecourse , ti//dt_timecourse*dt_timecourse + dt_timecourse 
                    if len(list(df[df['t'] <= tmin]['t'])) == 0:
                        tmin_event = tmin
                    else:
                        tmin_event = list(df[df['t'] <= tmin]['t'])[-1]
                        
                    if len(list(df[df['t'] >= tmax]['t'])) == 0:
                        tmax_event = tmax
                    else:
                        tmax_event = list(df[df['t'] >= tmax]['t'])[0]
                        
                    df_show = df[(df['t'] >= tmin_event) & (df['t'] <= tmax_event)]
                """
                
                if ti < dt_timecourse/2:
                    tmin_event = 0
                    tmax_event = tmax           
                else:
                    tmin_event = ti - dt_timecourse
                    tmax_event = ti + dt_timecourse  
                df_show = df[(df['t'] >= tmin_event) & (df['t'] <= tmax_event)]
                
                
                # v03 eventをすぎた段階でインデックスを更新する
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

                #color bar設定
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
                
                """
                # 真っ白なカラーバーを設定し、AX1と並ぶようにする
                for AX in [AX2,AX3]:
                    AX_inset = inset_axes(AX,
                                        width="100%",  # width = 50% of parent_bbox width
                                        height="5%",  # height : 5%
                                        bbox_to_anchor=(-0.025, 1.05, 1, 1),
                                        bbox_transform=AX.transAxes,
                                        loc='lower left')
                    cbar = fig.colorbar(cm.ScalarMappable(norm=None, cmap=ListedColormap(['w','w'])), cax=AX_inset,shrink=0.7, orientation="horizontal", drawedges=True)
                    cbar.outline.set_color('w')            
                    cbar.dividers.set_color('w')
                    cbar.ax.axis('off')
                """

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


            '''
            
            for i in range(len(df)):
                t  = df['t'].iloc[i]
                cx = df['X (pixel)'].iloc[i] 
                cy = df['Y (pixel)'].iloc[i] 
                cx_list =  df['X (pixel)'].iloc[:i+1]
                cy_list =  df['Y (pixel)'].iloc[:i+1] 
    
                event_position = tuple([int(df['ex'].iloc[i]), int(df['ey'].iloc[i])])
                if df['event'].iloc[i] == 1:                
                    RNA_state[event_position] = 1
                if df['event'].iloc[i] == 3:                
                    RNA_state[event_position] = 2
    
    
                xmin_disp, xmax_disp= min([cx-display_size,min(cx_list) ]), max([cx+display_size,max(cx_list) ])
                ymin_disp, ymax_disp= min([cy-display_size,min(cy_list) ]), max([cy+display_size,max(cy_list) ])
                min_disp, max_disp = min([xmin_disp,ymin_disp, min_disp]), max([xmax_disp,ymax_disp, max_disp])
                min_disp_nm, max_disp_nm = (min_disp-RNA_size//2)*pixel_to_nm, (max_disp-RNA_size//2)*pixel_to_nm 


                fig, ax = plt.subplots(figsize=(17, 10))
                gs = GridSpec(3, 3, left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.35, hspace=0.05)
                AX1 = plt.subplot(gs[:2, 0])
                AX2 = plt.subplot(gs[:2, 1])
                AX3 = plt.subplot(gs[:2, 2])
                AX4 = plt.subplot(gs[2:, :])
            
            
                for AXi in [AX2, AX3]:                      
                    AXi.text(min_disp_nm+(max_disp_nm-min_disp_nm)*0.65 ,max_disp_nm-(max_disp_nm-min_disp_nm)*0.85, 'Event = '+str(i)+'\n'+'t = '+ str(round_sig(t,3)) + ' s' ,color = 'k')
            
                cx_nm = (cx-RNA_size/2)*pixel_to_nm 
                cy_nm = (cy-RNA_size/2)*pixel_to_nm 
                cx_nm_list = (np.array(cx_list)-RNA_size/2) * pixel_to_nm 
                cy_nm_list = (np.array(cy_list)-RNA_size/2) * pixel_to_nm 

                #######################################################################################################################
                #######################################################################################################################
           #     AX1.set_title("p" + pi + ', [RNaseH]='+ str(RNaseH) + 'nM'+'\n')
                zoom_size = 20
                zoom_size_nm = zoom_size * pixel_to_nm
                AX1.text(cx_nm+zoom_size_nm*0.30 ,cy_nm-zoom_size_nm*0.7, 'Event = '+str(i)+'\n'+'t = '+ str(round_sig(t,3)) + ' s' ,color = 'k', ha='left')
                scalebar = 20 # nm
                AX1.plot( [cx_nm+zoom_size_nm*0.5,cx_nm+zoom_size_nm*0.5 + scalebar],[cy_nm+zoom_size_nm*0.7,cy_nm+zoom_size_nm*0.7], lw = 8,c='k')
                AX1.text(cx_nm+zoom_size_nm*0.5+ scalebar/2 ,cy_nm+zoom_size_nm*0.83, str(scalebar)+' nm', va='center', ha='center',weight='bold',fontsize= 20,zorder=10)
                m1 = AX1.imshow(RNA_state[cx-zoom_size:cx+zoom_size,cy-zoom_size:cy+zoom_size].T, clim=(0,4),cmap=colortable_RNA\
                           ,extent = ( cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm, cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm), origin='lower')    
                AX1.set(xlim=(cx_nm-zoom_size*pixel_to_nm, cx_nm+zoom_size*pixel_to_nm),ylim=(cy_nm-zoom_size*pixel_to_nm, cy_nm+zoom_size*pixel_to_nm) )

                #color bar設定
                AX1_inset = inset_axes(AX1,
                                    width="100%",  # width = 50% of parent_bbox width
                                    height="5%",  # height : 5%
                                    bbox_to_anchor=(-0.025, 1.05, 1, 1),
                                    bbox_transform=AX1.transAxes,
                                    loc='lower left')
                cbar = fig.colorbar(m1, cax=AX1_inset, orientation="horizontal",ticks=[0.5,1.5,2.5,3.5])
                cbar.ax.set_xticklabels(['Unbound', 'Bound', 'Hydrolyzed', 'Empty'])
                AX1_inset.xaxis.set_ticks_position("top")
                #######################################################################################################################
                #######################################################################################################################

                m2 = AX2.imshow(RNA_state[min_disp:max_disp,min_disp:max_disp].T, clim=(0,4),cmap=colortable_RNA\
                           ,extent = ( min_disp_nm, max_disp_nm, min_disp_nm, max_disp_nm), origin='lower')    
                

                # 真っ白なカラーバーを設定し、AX1と並ぶようにする
                for AX in [AX2,AX3]:
                    AX_inset = inset_axes(AX,
                                        width="100%",  # width = 50% of parent_bbox width
                                        height="5%",  # height : 5%
                                        bbox_to_anchor=(-0.025, 1.05, 1, 1),
                                        bbox_transform=AX.transAxes,
                                        loc='lower left')
                    cbar = fig.colorbar(cm.ScalarMappable(norm=None, cmap=ListedColormap(['w','w'])), cax=AX_inset,shrink=0.7, orientation="horizontal", drawedges=True)
                    cbar.outline.set_color('w')            
                    cbar.dividers.set_color('w')
                    cbar.ax.axis('off')


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
    
                imin,imax =( i//1000)*1000, min( [len(df)-1,(i//1000+1)*1000]) 
                
                AX4.axvline(x=df['t'].iloc[i], ls = 'dashed',c='0.5', lw = 1)
                
                AX4.plot([df['t'].iloc[i]], [74], marker= "v",c='0.5',ms=12, clip_on=False)
                AX4.plot([df['t'].iloc[i]], [df['N1'].iloc[i]], marker= "o",c='0.2',ms=10)
                AX4.plot([df['t'].iloc[i]], [ df['N2'].iloc[i] + df['N3'].iloc[i]], marker= "o",c='g',ms=10)
                AX4.plot([df['t'].iloc[i]], [df['N4'].iloc[i]], marker= "o",c='r',ms=10)
                
                AX4.step(df['t'].iloc[imin:imax],df['N1'].iloc[imin:imax],where='post',color = '0.2',label ='Unbound',lw=2)
                AX4.step(df['t'].iloc[imin:imax],df['N2'].iloc[imin:imax] + df['N3'].iloc[imin:imax],where='post',color = 'g',label ='Bound',lw=2)
                AX4.step(df['t'].iloc[imin:imax],df['N4'].iloc[imin:imax],where='post',color = 'r',label ='Hydrolyzed',lw=2)

                
                AX4.set(ylabel='Number of sites',xlabel='Time (s)',ylim=(0,70),yticks=[0,10,20,30,40,50,60,70],xlim=(df['t'].iloc[imin],df['t'].iloc[imax]) )            
                
                
                AX4.grid(which = 'both',axis='y')
                AX4.legend(loc ='upper right', ncol= 3)            
                
                
                #plt.subplots_adjust(hspace=0.05, wspace=0.05) 
                plt.savefig(trace_fol + os.sep + pi +'_'+ ('00000'+str(i))[-6:] +'.png' ,dpi = 100)
                plt.show()
                if i%10 == 0:
                    print('############################################################')
                    print(data)
                    print ('i='+str(i), ' t='+str(t))
                    plt.show()
                    MEMORY_RELEASE()
                plt.clf()
                plt.close()


            '''





    
if MAKE_MOVIE:
    for path in path_list:
        datas = sorted(glob.glob(path + os.sep + 'progress' + os.sep + '*.csv'))
        datas =  [f for f in datas if f[-8:-4].isnumeric()]
    
    
        for data in datas:
            pi = data[-8:-4]
            
            #command = r'ffmpeg -y -framerate 100 -i "'+trace_i_str+os.sep+trace_i_str+ '_%06d.png" -c:v libx264 -preset slow -crf 22 '+trace_i_str+ '.avi'
          #  command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -vf crop=iw:680:0:0 -pix_fmt yuv420p -c:v libx264 -preset slow  -crf 18 '+pi+ '.mp4'            
           # command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -c:v libx265 -preset slow  -crf 18 '+pi+ '.mp4'            
            command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -q 8 -crf 0 '+pi+ '.avi' 
            command = r'ffmpeg -y -framerate 100 -i "'+pi+os.sep+pi+ '_%06d.png" -pix_fmt yuv420p '+pi+ '.mp4' 
            print(command)
            p = subprocess.run(command, stdout=subprocess.PIPE, shell=True, cwd=path + os.sep + 'progress',universal_newlines=True)
        
        
    
    
    
    
    













