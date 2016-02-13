import matplotlib.pyplot as pl
from astropy.io import fits
from math import sqrt,ceil
from matplotlib import cm
import numpy as np
import glob as g
import pymysql
import os

pl.rcParams['legend.fontsize'] = 12

# to do:
# 	argparse for:
#		camera_id
#		dusty_centre
#		clean_centre
# 		box_width
#

camera_id=810
data_dir='/home/jmcc/data/ngts/paranal/dust/%d' % (camera_id)

def getData(camera_id):
	db=pymysql.connect(host='ngtsdb',db='ngts_ops')
	qry="SELECT action_id FROM action_summary_log WHERE camera_id=%d AND num_images>10 AND action_type='flatField'" % (camera_id)
	with db.cursor() as cur:
		cur.execute(qry)
		for row in cur:
			action_id=row[0]
			outdir='%s/%d/action%s_flatField/' % (data_dir,camera_id,action_id)
			if os.path.exists(outdir) == False:
				os.mkdir(outdir)
			qry2="SELECT image_id FROM raw_image_list WHERE action_id=%s AND exptime > 2" % (action_id)
			with db.cursor() as cur2:
				cur2.execute(qry2)
				for row2 in cur2:
					image_id=row2[0]
					comm='cp /ngts/testdata/paranal/action%s_flatField/IMAGE%s.fits.bz2 %s' % (action_id,image_id,outdir)
					print comm
					os.system(comm)
	db.close()

def makeMasterFlat(action):
	os.chdir(action)
	t2=g.glob('*.fits')
	dc_list=[]
	flat_combined_total=0
	for i in t2:
		h=fits.open(i)
		d=h[0].data[0:2048,20:2068]
		oscan=h[0].data[4:,-15:]
		dc=d-np.median(oscan)
		dc_list.append(dc)
		flat_combined_total+=1
		h.close()
	dc_cube=np.dstack(dc_list)
	masterFlat=np.median(dc_cube,axis=2)
	print masterFlat.shape, flat_combined_total
	os.chdir('../')
	return masterFlat

def compareDustyRegion(img):
	# coordinate assume the oscan has been removed
	width=250
	half_width=width/2
	dusty_centre_x=1114  # 1134 with oscan
	dusty_centre_y=639   #  659 with oscan
	clean_centre_x=780   #  800 with oscan
	clean_centre_y=1173  # 1193 with oscan
	clean_region=img[clean_centre_y-half_width:clean_centre_y+half_width,clean_centre_x-half_width:clean_centre_x+half_width]
	dusty_region=img[dusty_centre_y-half_width:dusty_centre_y+half_width,dusty_centre_x-half_width:dusty_centre_x+half_width]
	fr=np.median(dusty_region)/np.median(clean_region)
	return fr

def main():
	os.chdir(data_dir)
	action_list=g.glob('action*')
	print('Making master reference flat...')
	fr,fr_f=[],[]
	masterRefFlat=makeMasterFlat(action_list[0])
	fr.append(compareDustyRegion(masterRefFlat))
	
	# compare the ratio of flux in and out of the dusty region
	# and make some plots of the flat fielded master flats too
	fig = pl.figure(1,figsize=(15,15))
	montage_dim=int(ceil(sqrt(len(action_list)-1)))
	for i in range(1,len(action_list)):
		print('[%d/%d] Making master flat %s' % (i,(len(action_list)-1),action_list[i]))
		masterFlat=makeMasterFlat(action_list[i])
		ratio=masterFlat/masterRefFlat
		fr.append(compareDustyRegion(masterFlat))
		fr_f.append(compareDustyRegion(ratio))
		ax = fig.add_subplot(montage_dim, montage_dim, i, xticks=[], yticks=[])
		ax.imshow(ratio,cmap=cm.afmhot,vmin=0.8*np.median(ratio),vmax=1.2*np.median(ratio),interpolation=None)        
		ax.set_title('%s' % (action_list[i].split('_')[0].split('action')[1]))
	print('Saving figures...')
	fig.savefig('dustEvolution_imgs.png',dpi=200)	

	fig2,ax2=pl.subplots(2,1,figsize=(10,10))
	ax2[0].plot(fr,'ro')
	ax2[0].set_xlabel('Time (action_id)')
	ax2[0].set_ylabel('dusty/clean flux ratio')
	ax2[1].plot(fr_f,'bo')
	ax2[1].set_xlabel('Time (action_id)')
	ax2[1].set_ylabel('dusty/clean flux ratio (ff)')
	fig2.savefig('dustEvolution_graph.png',dpi=200)

if __name__ == '__main__':
	main()

