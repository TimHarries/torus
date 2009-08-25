#!/bin/ksh
#set -ax

#################################################
#  Grid generation study
#################################################

############
# Figure 1 #
############

# Get mesh plots
scp zen:/scratch/acreman/sph_disc/grid_tests/working/runs/cyl_geom_F/npart_1e7/set1/run_5.e26/mesh0000.tif .
/sw/bin/convert mesh0000.tif set1_mesh.pdf
rm mesh0000.tif

scp zen:/scratch/acreman/sph_disc/grid_tests/working/runs/cyl_geom_F/npart_1e7/set2/run_0.1/mesh0000.tif .
/sw/bin/convert mesh0000.tif set2_mesh.pdf
rm mesh0000.tif

############
# Figure 2 #
############

# Get plots from grid generation parameter study
# Cartesian grids
export BASE_DIR=zen:/scratch/acreman/sph_disc/grid_tests/working/runs/cyl_geom_F/
for num_part in 1e5 1e6  1e7; do
    work_dir=npart_${num_part}
    for dir in set1 set2; do
	for file in ${dir}.ps; do
	    this_fig=${dir}.ps
	    scp ${BASE_DIR}/${work_dir}/${dir}/${this_fig} ${num_part}_${this_fig}
	    ps2pdf ${num_part}_${this_fig}
	    rm ${num_part}_${this_fig}
	done
    done
done

############
# Figure 3 #
############

export BASE_DIR=zen:/scratch/acreman/sph_disc/grid_tests/working/runs/cyl_geom_T/
for num_part in 1e5 1e6  1e7; do
    work_dir=npart_${num_part}
    for dir in set1 set2; do
	for file in ${dir}.ps; do
	    this_fig=${dir}.ps
	    scp ${BASE_DIR}/${work_dir}/${dir}/${this_fig} ${num_part}_cyl_${this_fig}
	    ps2pdf ${num_part}_cyl_${this_fig}
	    rm ${num_part}_cyl_${this_fig}
	done
    done
done

#################################################
# sphbench
#################################################

############
# Figure 4 #
############

# Get SEDs for 10^5 particles and no grid modifications. 
scp zen:/scratch/acreman/sph_disc/sphbench/no_mod_runs/no_mod/sed_013.ps .
scp zen:/scratch/acreman/sph_disc/sphbench/no_mod_runs/no_mod/sed_077.ps .
ps2pdf sed_013.ps
ps2pdf sed_077.ps
rm sed_013.ps sed_077.ps

############
# Figure 5 #
############

scp zen:/scratch/acreman/sph_disc/sphbench/no_mod_runs/no_mod/visit0000.tif .
/sw/bin/convert visit0000.tif basic_dens.pdf
rm  visit0000.tif

############
# Figure 6 #
############

scp zen:/scratch/acreman/sph_disc/sphbench/forced/sed_013.ps  seds_forced_13.ps
scp zen:/scratch/acreman/sph_disc/sphbench/forced/sed_077.ps  seds_forced_77.ps
ps2pdf seds_forced_13.ps
ps2pdf seds_forced_77.ps
rm seds_forced_13.ps seds_forced_77.ps

############
# Figure 7 #
############

# Density plots for different numbers of particles
scp zen:/scratch/acreman/sph_disc/sphbench/npart/1e5/visit0000.tif .
/sw/bin/convert visit0000.tif dens_1e5.pdf
rm visit0000.tif

scp zen:/scratch/acreman/sph_disc/sphbench/npart/1e6/visit0000.tif .
/sw/bin/convert visit0000.tif dens_1e6.pdf
rm visit0000.tif

scp zen:/scratch/acreman/sph_disc/sphbench/npart/1e7/visit0000.tif .
/sw/bin/convert visit0000.tif dens_1e7.pdf
rm visit0000.tif

scp zen:/scratch/acreman/sph_disc/sphbench/npart/1e8/visit0000.tif .
/sw/bin/convert visit0000.tif dens_1e8.pdf
rm visit0000.tif

############
# Figure 8 #
############

# Temperature differences
for thisRun in 1e5 1e6 1e7 1e8; do
    scp zen:/scratch/acreman/sph_disc/sphbench/npart/${thisRun}_tdiff/pgplot.ps tdiff_${thisRun}.ps
    ps2pdf tdiff_${thisRun}.ps
done

############
# Figure 9 #
############

# Midplane temperature differences
for  thisRun in 1 2 3 4; do
    scp zen:/scratch/acreman/sph_disc/sphbench/npart/mp_comp/run/mp_temp_${thisRun}.ps .
    ps2pdf mp_temp_${thisRun}.ps
    rm mp_temp_${thisRun}.ps
done

###################
# Figures 10      #
###################

# SEDs for different numbers of particles. 
for run in 1e5 1e6 1e7 1e8; do
    scp zen:/scratch/acreman/sph_disc/sphbench/npart/${run}/sed_013.ps sed_013_${run}.ps
    scp zen:/scratch/acreman/sph_disc/sphbench/npart/${run}/sed_077.ps sed_077_${run}.ps
done

for sed_plot in sed_0??_1e?.ps
do
    ps2pdf ${sed_plot}
    rm ${sed_plot}
done

#################################################
#  SPH-Torus disc
#################################################

# Plots showing adjustment to equilibrium
scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run/dumps/large/u_DSC6675.ps . 
ps2pdf u_DSC6675.ps
rm u_DSC6675.ps

scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run/dumps/large/u_DSC6691.ps . 
ps2pdf u_DSC6691.ps
rm u_DSC6691.ps

scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run/dumps/large/u_DSC6763.ps . 
ps2pdf u_DSC6763.ps
rm u_DSC6763.ps

scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run/dumps/large/u_DSC7363.ps .
ps2pdf u_DSC7363.ps
rm u_DSC7363.ps

scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run/dumps/large/u_DSC9363.ps .
ps2pdf u_DSC9363.ps
rm u_DSC9363.ps

scp zen:/scratch/acreman/SPH-Torus_runs/run0001/run2/dumps/u_run2_DSC1024.ps .
ps2pdf u_run2_DSC1024.ps
rm u_run2_DSC1024.ps

cp ~/sph_disc_pp/profiles/sigma.ps .
ps2pdf sigma.ps

cp ~/sph_disc_pp/profiles/sigma_final.ps .
ps2pdf sigma_final.ps

cp ~/sph_disc_pp/scale_height/mid_tem.ps .
ps2pdf mid_tem.ps

cp ~/sph_disc_pp/scale_height/scale_height.ps .
ps2pdf scale_height.ps

rm *.ps

exit

