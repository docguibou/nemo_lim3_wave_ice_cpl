# Automatic Makefile

FCM_BLD_CFG = /home1/datawork/gboutin/trunk/NEMOGCM/TOOLS/COMPILE/bldxag.cfg
export FCM_VERBOSE ?= 1

export FCM_ROOTDIR := /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/BLD
export FCM_CACHEDIR := $(FCM_ROOTDIR)/.cache/.bld
export FCM_CFGDIR := $(FCM_ROOTDIR)/cfg
export FCM_SRCDIR := $(FCM_ROOTDIR)/src
export FCM_BINDIR := $(FCM_ROOTDIR)/bin
export FCM_ETCDIR := $(FCM_ROOTDIR)/etc
export FCM_LIBDIR := $(FCM_ROOTDIR)/lib
export FCM_TMPDIR := $(FCM_ROOTDIR)/tmp
export FCM_DONEDIR := $(FCM_ROOTDIR)/done
export FCM_FLAGSDIR := $(FCM_ROOTDIR)/flags
export FCM_INCDIR := $(FCM_ROOTDIR)/inc
export FCM_PPSRCDIR := $(FCM_ROOTDIR)/ppsrc
export FCM_OBJDIR := $(FCM_ROOTDIR)/obj
export FCM_ROOTPATH := $(FCM_ROOTDIR)
export FCM_CACHEPATH := $(FCM_CACHEDIR)
export FCM_CFGPATH := $(FCM_CFGDIR)
export FCM_SRCPATH := $(FCM_SRCDIR)
export FCM_BINPATH := $(FCM_BINDIR)
export FCM_ETCPATH := $(FCM_ETCDIR)
export FCM_LIBPATH := $(FCM_LIBDIR)
export FCM_TMPPATH := $(FCM_TMPDIR)
export FCM_DONEPATH := $(FCM_DONEDIR)
export FCM_FLAGSPATH := $(FCM_FLAGSDIR)
export FCM_INCPATH := $(FCM_INCDIR)
export FCM_PPSRCPATH := $(FCM_PPSRCDIR)
export FCM_OBJPATH := $(FCM_OBJDIR)

export PERL5LIB := /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/fcm/bin/../lib

export AR = ar
export ARFLAGS = -r
export CC = cc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = 
export CPPKEYS = 
export DIFF3FLAGS = -E -m
export FC = mpif90
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_MODSEARCH = 
export FC_OUTPUT = -o
export FFLAGS = -O2 -fno-alias -fno-fnalias -r8 -ip -ftz -sox -assume byterecl -convert big_endian -I/home1/datahome/gboutin/XIOS/xios-2.0-mpt/inc -I/appli/NETCDF/4.3.3.1-mpt-intel/include -I/home1/datahome/gboutin/OASIS3-MCT_2.0_branch/oasis3-mct/compile_oa3-mct_gfopenmpi/build/lib/psmile.MPI1 -I/home1/datahome/gboutin/OASIS3-MCT_2.0_branch/oasis3-mct/compile_oa3-mct_gfopenmpi/build/lib/mct
export FPPKEYS = key_lim3 key_mpp_mpi key_nosignedzero key_iomput key_oasis3
export LD = mpif90 -lstdc++ -lmpi
export LDFLAGS =  -L/home1/datahome/gboutin/XIOS/xios-2.0-mpt/lib -lxios -L/appli/NETCDF/4.3.3.1-mpt-intel/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lhdf5 -lz -lcurl -L/home1/datahome/gboutin/OASIS3-MCT_2.0_branch/oasis3-mct/compile_oa3-mct_gfopenmpi/lib -lpsmile.MPI1 -lmct -lmpeu -lscrip
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/sh

vpath % $(FCM_BINPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.etc $(FCM_ETCPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath % $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)

FCM_BLD_TARGETS = nemo.exe

.PHONY : all

all : $(FCM_BLD_TARGETS)

$(FCM_DONEDIR)/FCM_CP.dummy:
	touch $@

CFLAGS__nemo.flags: CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__nemo.flags: CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__nemo.flags: FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo.flags: FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo.flags: LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo.flags: LD.flags
	touch $(FCM_FLAGSDIR)/$@

OBJECTS__nemo = albedo.o asmbkg.o asminc.o asmpar.o bdy_oce.o bdydta.o bdydyn.o bdydyn2d.o bdydyn3d.o bdyice_lim.o bdyini.o bdylib.o bdytides.o bdytra.o bdyvol.o c1d.o cool_skin.o cpl_oasis3.o crs.o crsdom.o crsdomwri.o crsfld.o crsini.o crslbclnk.o cyclone.o daymod.o depth_e3.o dia25h.o diaar5.o diacfl.o diadct.o diaharm.o diahsb.o diahth.o dianam.o diaobs.o diaptr.o diatmb.o diawri.o diurnal_bulk.o divhor.o dom_oce.o domain.o domc1d.o domhgr.o dommsk.o domngb.o domvvl.o domwri.o domzgr.o dtatsd.o dtauvd.o dynadv.o dynadv_cen2.o dynadv_ubs.o dynbfr.o dyncor_c1d.o dyndmp.o dynhpg.o dynkeg.o dynldf.o dynldf_iso.o dynldf_lap_blp.o dynnxt.o dynspg.o dynspg_exp.o dynspg_ts.o dynvor.o dynzad.o dynzdf.o dynzdf_exp.o dynzdf_imp.o eosbn2.o fldread.o flo4rk.o flo_oce.o floats.o floblk.o flodom.o florst.o flowri.o geo2ocean.o icb_oce.o icbclv.o icbdia.o icbdyn.o icbini.o icblbc.o icbrst.o icbstp.o icbthm.o icbtrj.o icbutl.o ice.o in_out_manager.o iom.o iom_def.o iom_nf90.o iscplhsb.o iscplini.o iscplrst.o istate.o julian.o lbclnk.o lbcnfd.o ldfc1d_c2d.o ldfdyn.o ldfslp.o ldftra.o lib_cray.o lib_fortran.o lib_mpp.o limadv_prather.o limadv_umx.o limcons.o limctl.o limdiahsb.o limdyn.o limhdf.o limistate.o limitd_me.o limitd_th.o limrhg.o limrst.o limsbc.o limtab.o limthd.o limthd_da.o limthd_dh.o limthd_dif.o limthd_ent.o limthd_lac.o limthd_sal.o limtrp.o limupdate1.o limupdate2.o limvar.o limwri.o mpp_map.o mppini.o nemogcm.o obs_const.o obs_conv.o obs_fbm.o obs_grid.o obs_inter_h2d.o obs_inter_sup.o obs_inter_z1d.o obs_mpp.o obs_oper.o obs_prep.o obs_profiles.o obs_profiles_def.o obs_read_altbias.o obs_read_prof.o obs_read_surf.o obs_readmdt.o obs_rot_vel.o obs_sort.o obs_sstbias.o obs_surf_def.o obs_types.o obs_utils.o obs_write.o oce.o par_kind.o par_oce.o phycst.o prtctl.o restart.o sbc_ice.o sbc_oce.o sbcapr.o sbcblk.o sbcblk_algo_coare.o sbcblk_algo_coare3p5.o sbcblk_algo_ecmwf.o sbcblk_algo_ncar.o sbccpl.o sbcdcy.o sbcflx.o sbcfwb.o sbcice_cice.o sbcice_if.o sbcice_lim.o sbcice_lim_2.o sbcisf.o sbcmod.o sbcrnf.o sbcssm.o sbcssr.o sbctide.o sbcwave.o solfrac_mod.o sshwzv.o step.o step_c1d.o step_diu.o step_oce.o stopar.o stopts.o storng.o stpctl.o thd_ice.o tide_mod.o tideini.o timing.o traadv.o traadv_cen.o traadv_fct.o traadv_mle.o traadv_mus.o traadv_qck.o traadv_ubs.o trabbc.o trabbl.o tradmp.o traldf.o traldf_iso.o traldf_lap_blp.o traldf_triad.o tranpc.o tranxt.o traqsr.o trasbc.o trazdf.o trazdf_exp.o trazdf_imp.o trc_oce.o trd_oce.o trddyn.o trdglo.o trdini.o trdken.o trdmxl.o trdmxl_oce.o trdmxl_rst.o trdpen.o trdtra.o trdtrc.o trdvor.o trdvor_oce.o updtide.o usrdef_closea.o usrdef_fmask.o usrdef_hgr.o usrdef_istate.o usrdef_nam.o usrdef_sbc.o usrdef_zgr.o wet_dry.o wrk_nemo.o zdf_oce.o zdfbfr.o zdfddm.o zdfevd.o zdfgls.o zdfini.o zdfmxl.o zdfqiao.o zdfric.o zdftke.o zdftmx.o zpshde.o

libnemo.a: $(OBJECTS__nemo)
	fcm_internal archive $@ $^

CFLAGS__ioipsl.flags: CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__ioipsl.flags: CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__ioipsl.flags: FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__ioipsl.flags: FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__ioipsl.flags: LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

LD__ioipsl.flags: LD.flags
	touch $(FCM_FLAGSDIR)/$@

OBJECTS__ioipsl = calendar.o defprec.o errioipsl.o flincom.o fliocom.o getincom.o histcom.o ioipsl.o mathelp.o nc4interface.o restcom.o stringop.o

libioipsl.a: $(OBJECTS__ioipsl)
	fcm_internal archive $@ $^

CC.flags:
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags: CC.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

FC.flags:
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags: FC.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

LD.flags:
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags:
	touch $(FCM_FLAGSDIR)/$@

export OBJECTS = $(OBJECTS__ioipsl) $(OBJECTS__nemo)

libfcm_default.a: $(OBJECTS)
	fcm_internal archive $@ $^

FFLAGS__ioipsl__calendar.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

calendar.done: calendar.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

calendar.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/calendar.f90 FFLAGS__ioipsl__calendar.flags errioipsl.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__defprec.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

defprec.done: defprec.o
	touch $(FCM_DONEDIR)/$@

defprec.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/defprec.f90 FFLAGS__ioipsl__defprec.flags
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__errioipsl.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

errioipsl.done: errioipsl.o
	touch $(FCM_DONEDIR)/$@

errioipsl.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/errioipsl.f90 FFLAGS__ioipsl__errioipsl.flags
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__flincom.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

flincom.done: flincom.o calendar.done errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

flincom.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/flincom.f90 FFLAGS__ioipsl__flincom.flags calendar.o errioipsl.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__fliocom.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

fliocom.done: fliocom.o calendar.done defprec.done errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

fliocom.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/fliocom.f90 FFLAGS__ioipsl__fliocom.flags calendar.o defprec.o errioipsl.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__getincom.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

getincom.done: getincom.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

getincom.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/getincom.f90 FFLAGS__ioipsl__getincom.flags errioipsl.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__histcom.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

histcom.done: histcom.o calendar.done errioipsl.done fliocom.done mathelp.done nc4interface.done stringop.done
	touch $(FCM_DONEDIR)/$@

histcom.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/histcom.f90 FFLAGS__ioipsl__histcom.flags calendar.o errioipsl.o fliocom.o mathelp.o nc4interface.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__ioipsl.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl.done: ioipsl.o calendar.done errioipsl.done flincom.done fliocom.done getincom.done histcom.done mathelp.done restcom.done stringop.done
	touch $(FCM_DONEDIR)/$@

ioipsl.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/ioipsl.f90 FFLAGS__ioipsl__ioipsl.flags calendar.o errioipsl.o flincom.o fliocom.o getincom.o histcom.o mathelp.o restcom.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__mathelp.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

mathelp.done: mathelp.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

mathelp.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/mathelp.f90 FFLAGS__ioipsl__mathelp.flags errioipsl.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__nc4interface.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

nc4interface.done: nc4interface.o
	touch $(FCM_DONEDIR)/$@

nc4interface.o: $(FCM_PPSRCDIR)/ioipsl/nc4interface.f90 FFLAGS__ioipsl__nc4interface.flags
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__restcom.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

restcom.done: restcom.o calendar.done errioipsl.done fliocom.done mathelp.done stringop.done
	touch $(FCM_DONEDIR)/$@

restcom.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/restcom.f90 FFLAGS__ioipsl__restcom.flags calendar.o errioipsl.o fliocom.o mathelp.o stringop.o
	fcm_internal compile:F ioipsl $< $@

FFLAGS__ioipsl__stringop.flags: FFLAGS__ioipsl.flags
	touch $(FCM_FLAGSDIR)/$@

stringop.done: stringop.o
	touch $(FCM_DONEDIR)/$@

stringop.o: /home1/datawork/gboutin/trunk/NEMOGCM/EXTERNAL/IOIPSL/src/stringop.f90 FFLAGS__ioipsl__stringop.flags
	fcm_internal compile:F ioipsl $< $@

FFLAGS__nemo__albedo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

albedo.done: albedo.o in_out_manager.done lib_fortran.done lib_mpp.done phycst.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

albedo.o: $(FCM_PPSRCDIR)/nemo/albedo.f90 FFLAGS__nemo__albedo.flags in_out_manager.o lib_fortran.o lib_mpp.o phycst.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__asmbkg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

asmbkg.done: asmbkg.o asmpar.done dom_oce.done eosbn2.done ice.done in_out_manager.done iom.done ldfslp.done ldftra.done oce.done sbc_oce.done tradmp.done zdf_oce.done zdfddm.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

asmbkg.o: $(FCM_PPSRCDIR)/nemo/asmbkg.f90 FFLAGS__nemo__asmbkg.flags asmpar.o dom_oce.o eosbn2.o ice.o in_out_manager.o iom.o ldfslp.o ldftra.o oce.o sbc_oce.o tradmp.o zdf_oce.o zdfddm.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__asminc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

asminc.done: asminc.o asmpar.done c1d.done diaobs.done dom_oce.done domvvl.done eosbn2.done in_out_manager.done iom.done ldfdyn.done lib_mpp.done oce.done par_oce.done sbc_oce.done wrk_nemo.done zpshde.done
	touch $(FCM_DONEDIR)/$@

asminc.o: $(FCM_PPSRCDIR)/nemo/asminc.f90 FFLAGS__nemo__asminc.flags asmpar.o c1d.o diaobs.o dom_oce.o domvvl.o eosbn2.o in_out_manager.o iom.o ldfdyn.o lib_mpp.o oce.o par_oce.o sbc_oce.o wrk_nemo.o zpshde.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__asmpar.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

asmpar.done: asmpar.o
	touch $(FCM_DONEDIR)/$@

asmpar.o: $(FCM_PPSRCDIR)/nemo/asmpar.f90 FFLAGS__nemo__asmpar.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdy_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdy_oce.done: bdy_oce.o lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

bdy_oce.o: $(FCM_PPSRCDIR)/nemo/bdy_oce.f90 FFLAGS__nemo__bdy_oce.flags lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdydta.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdydta.done: bdydta.o bdy_oce.done bdytides.done dom_oce.done fldread.done ice.done in_out_manager.done iom.done limvar.done oce.done phycst.done sbcapr.done sbctide.done timing.done
	touch $(FCM_DONEDIR)/$@

bdydta.o: $(FCM_PPSRCDIR)/nemo/bdydta.f90 FFLAGS__nemo__bdydta.flags bdy_oce.o bdytides.o dom_oce.o fldread.o ice.o in_out_manager.o iom.o limvar.o oce.o phycst.o sbcapr.o sbctide.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdydyn.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdydyn.done: bdydyn.o bdy_oce.done bdydyn2d.done bdydyn3d.done dom_oce.done domvvl.done in_out_manager.done lbclnk.done oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

bdydyn.o: $(FCM_PPSRCDIR)/nemo/bdydyn.f90 FFLAGS__nemo__bdydyn.flags bdy_oce.o bdydyn2d.o bdydyn3d.o dom_oce.o domvvl.o in_out_manager.o lbclnk.o oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdydyn2d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdydyn2d.done: bdydyn2d.o bdy_oce.done bdylib.done dom_oce.done in_out_manager.done lbclnk.done oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

bdydyn2d.o: $(FCM_PPSRCDIR)/nemo/bdydyn2d.f90 FFLAGS__nemo__bdydyn2d.flags bdy_oce.o bdylib.o dom_oce.o in_out_manager.o lbclnk.o oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdydyn3d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdydyn3d.done: bdydyn3d.o bdy_oce.done bdylib.done dom_oce.done in_out_manager.done lbclnk.done oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

bdydyn3d.o: $(FCM_PPSRCDIR)/nemo/bdydyn3d.f90 FFLAGS__nemo__bdydyn3d.flags bdy_oce.o bdylib.o dom_oce.o in_out_manager.o lbclnk.o oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdyice_lim.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdyice_lim.done: bdyice_lim.o bdy_oce.done dom_oce.done eosbn2.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limctl.done limvar.done oce.done par_oce.done phycst.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

bdyice_lim.o: $(FCM_PPSRCDIR)/nemo/bdyice_lim.f90 FFLAGS__nemo__bdyice_lim.flags bdy_oce.o dom_oce.o eosbn2.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limctl.o limvar.o oce.o par_oce.o phycst.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdyini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdyini.done: bdyini.o bdy_oce.done bdydta.done bdytides.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done sbctide.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

bdyini.o: $(FCM_PPSRCDIR)/nemo/bdyini.f90 FFLAGS__nemo__bdyini.flags bdy_oce.o bdydta.o bdytides.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o sbctide.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdylib.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdylib.done: bdylib.o bdy_oce.done dom_oce.done in_out_manager.done lbclnk.done oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

bdylib.o: $(FCM_PPSRCDIR)/nemo/bdylib.f90 FFLAGS__nemo__bdylib.flags bdy_oce.o dom_oce.o in_out_manager.o lbclnk.o oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdytides.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdytides.done: bdytides.o bdy_oce.done daymod.done dom_oce.done fldread.done in_out_manager.done iom.done lbclnk.done oce.done phycst.done tideini.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

bdytides.o: $(FCM_PPSRCDIR)/nemo/bdytides.f90 FFLAGS__nemo__bdytides.flags bdy_oce.o daymod.o dom_oce.o fldread.o in_out_manager.o iom.o lbclnk.o oce.o phycst.o tideini.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdytra.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdytra.done: bdytra.o bdy_oce.done bdylib.done dom_oce.done in_out_manager.done lbclnk.done oce.done timing.done
	touch $(FCM_DONEDIR)/$@

bdytra.o: $(FCM_PPSRCDIR)/nemo/bdytra.f90 FFLAGS__nemo__bdytra.flags bdy_oce.o bdylib.o dom_oce.o in_out_manager.o lbclnk.o oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__bdyvol.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

bdyvol.done: bdyvol.o bdy_oce.done dom_oce.done in_out_manager.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done sbcisf.done timing.done
	touch $(FCM_DONEDIR)/$@

bdyvol.o: $(FCM_PPSRCDIR)/nemo/bdyvol.f90 FFLAGS__nemo__bdyvol.flags bdy_oce.o dom_oce.o in_out_manager.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o sbcisf.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__c1d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

c1d.done: c1d.o par_kind.done
	touch $(FCM_DONEDIR)/$@

c1d.o: $(FCM_PPSRCDIR)/nemo/c1d.f90 FFLAGS__nemo__c1d.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__cool_skin.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

cool_skin.done: cool_skin.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done par_kind.done phycst.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

cool_skin.o: $(FCM_PPSRCDIR)/nemo/cool_skin.f90 FFLAGS__nemo__cool_skin.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o par_kind.o phycst.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__cpl_oasis3.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

cpl_oasis3.done: cpl_oasis3.o dom_oce.done in_out_manager.done lbclnk.done par_oce.done
	touch $(FCM_DONEDIR)/$@

cpl_oasis3.o: $(FCM_PPSRCDIR)/nemo/cpl_oasis3.f90 FFLAGS__nemo__cpl_oasis3.flags dom_oce.o in_out_manager.o lbclnk.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crs.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crs.done: crs.o dom_oce.done in_out_manager.done par_oce.done
	touch $(FCM_DONEDIR)/$@

crs.o: $(FCM_PPSRCDIR)/nemo/crs.f90 FFLAGS__nemo__crs.flags dom_oce.o in_out_manager.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crsdom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crsdom.done: crsdom.o crs.done crslbclnk.done dom_oce.done in_out_manager.done lib_mpp.done par_kind.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

crsdom.o: $(FCM_PPSRCDIR)/nemo/crsdom.f90 FFLAGS__nemo__crsdom.flags crs.o crslbclnk.o dom_oce.o in_out_manager.o lib_mpp.o par_kind.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crsdomwri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crsdomwri.done: crsdomwri.o crs.done crsdom.done crslbclnk.done dom_oce.done in_out_manager.done iom.done iom_def.done lib_mpp.done par_kind.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

crsdomwri.o: $(FCM_PPSRCDIR)/nemo/crsdomwri.f90 FFLAGS__nemo__crsdomwri.flags crs.o crsdom.o crslbclnk.o dom_oce.o in_out_manager.o iom.o iom_def.o lib_mpp.o par_kind.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crsfld.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crsfld.done: crsfld.o crs.done crsdom.done crslbclnk.done dom_oce.done in_out_manager.done iom.done lbclnk.done ldftra.done oce.done sbc_oce.done timing.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

crsfld.o: $(FCM_PPSRCDIR)/nemo/crsfld.f90 FFLAGS__nemo__crsfld.flags crs.o crsdom.o crslbclnk.o dom_oce.o in_out_manager.o iom.o lbclnk.o ldftra.o oce.o sbc_oce.o timing.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crsini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crsini.done: crsini.o crs.done crsdom.done crsdomwri.done crslbclnk.done dom_oce.done in_out_manager.done iom.done lib_mpp.done par_kind.done par_oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

crsini.o: $(FCM_PPSRCDIR)/nemo/crsini.f90 FFLAGS__nemo__crsini.flags crs.o crsdom.o crsdomwri.o crslbclnk.o dom_oce.o in_out_manager.o iom.o lib_mpp.o par_kind.o par_oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__crslbclnk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

crslbclnk.done: crslbclnk.o crs.done dom_oce.done in_out_manager.done lbclnk.done par_kind.done
	touch $(FCM_DONEDIR)/$@

crslbclnk.o: $(FCM_PPSRCDIR)/nemo/crslbclnk.f90 FFLAGS__nemo__crslbclnk.flags crs.o dom_oce.o in_out_manager.o lbclnk.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__cyclone.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

cyclone.done: cyclone.o
	touch $(FCM_DONEDIR)/$@

cyclone.o: $(FCM_PPSRCDIR)/nemo/cyclone.f90 FFLAGS__nemo__cyclone.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__daymod.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

daymod.done: daymod.o dom_oce.done in_out_manager.done ioipsl.done iom.done phycst.done prtctl.done restart.done timing.done
	touch $(FCM_DONEDIR)/$@

daymod.o: $(FCM_PPSRCDIR)/nemo/daymod.f90 FFLAGS__nemo__daymod.flags dom_oce.o in_out_manager.o ioipsl.o iom.o phycst.o prtctl.o restart.o timing.o
	fcm_internal compile:F nemo $< $@

ddatetoymdhms.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/ddatetoymdhms.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

ddatetoymdhms.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/ddatetoymdhms.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__depth_e3.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

depth_e3.done: depth_e3.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

depth_e3.o: $(FCM_PPSRCDIR)/nemo/depth_e3.f90 FFLAGS__nemo__depth_e3.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dia25h.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dia25h.done: dia25h.o dom_oce.done in_out_manager.done iom.done oce.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

dia25h.o: $(FCM_PPSRCDIR)/nemo/dia25h.f90 FFLAGS__nemo__dia25h.flags dom_oce.o in_out_manager.o iom.o oce.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diaar5.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diaar5.done: diaar5.o dom_oce.done eosbn2.done fldread.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

diaar5.o: $(FCM_PPSRCDIR)/nemo/diaar5.f90 FFLAGS__nemo__diaar5.flags dom_oce.o eosbn2.o fldread.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diacfl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diacfl.done: diacfl.o dom_oce.done domvvl.done in_out_manager.done lbclnk.done lib_mpp.done oce.done timing.done
	touch $(FCM_DONEDIR)/$@

diacfl.o: $(FCM_PPSRCDIR)/nemo/diacfl.f90 FFLAGS__nemo__diacfl.flags dom_oce.o domvvl.o in_out_manager.o lbclnk.o lib_mpp.o oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diadct.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diadct.done: diadct.o
	touch $(FCM_DONEDIR)/$@

diadct.o: $(FCM_PPSRCDIR)/nemo/diadct.f90 FFLAGS__nemo__diadct.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diaharm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diaharm.done: diaharm.o
	touch $(FCM_DONEDIR)/$@

diaharm.o: $(FCM_PPSRCDIR)/nemo/diaharm.f90 FFLAGS__nemo__diaharm.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diahsb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diahsb.done: diahsb.o bdy_oce.done dom_oce.done domvvl.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done restart.done sbc_oce.done sbcisf.done sbcrnf.done timing.done trabbc.done traqsr.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

diahsb.o: $(FCM_PPSRCDIR)/nemo/diahsb.f90 FFLAGS__nemo__diahsb.flags bdy_oce.o dom_oce.o domvvl.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o restart.o sbc_oce.o sbcisf.o sbcrnf.o timing.o trabbc.o traqsr.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diahth.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diahth.done: diahth.o
	touch $(FCM_DONEDIR)/$@

diahth.o: $(FCM_PPSRCDIR)/nemo/diahth.f90 FFLAGS__nemo__diahth.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dianam.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dianam.done: dianam.o dom_oce.done in_out_manager.done ioipsl.done phycst.done
	touch $(FCM_DONEDIR)/$@

dianam.o: $(FCM_PPSRCDIR)/nemo/dianam.f90 FFLAGS__nemo__dianam.flags dom_oce.o in_out_manager.o ioipsl.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diaobs.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diaobs.done: diaobs.o dom_oce.done ice.done in_out_manager.done lib_mpp.done mpp_map.done obs_grid.done obs_oper.done obs_prep.done obs_profiles_def.done obs_read_altbias.done obs_read_prof.done obs_read_surf.done obs_readmdt.done obs_rot_vel.done obs_sstbias.done obs_surf_def.done obs_types.done obs_write.done oce.done par_kind.done par_oce.done phycst.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

diaobs.o: $(FCM_PPSRCDIR)/nemo/diaobs.f90 FFLAGS__nemo__diaobs.flags dom_oce.o ice.o in_out_manager.o lib_mpp.o mpp_map.o obs_grid.o obs_oper.o obs_prep.o obs_profiles_def.o obs_read_altbias.o obs_read_prof.o obs_read_surf.o obs_readmdt.o obs_rot_vel.o obs_sstbias.o obs_surf_def.o obs_types.o obs_write.o oce.o par_kind.o par_oce.o phycst.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diaptr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diaptr.done: diaptr.o dom_oce.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

diaptr.o: $(FCM_PPSRCDIR)/nemo/diaptr.f90 FFLAGS__nemo__diaptr.flags dom_oce.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diatmb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diatmb.done: diatmb.o dom_oce.done in_out_manager.done iom.done oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

diatmb.o: $(FCM_PPSRCDIR)/nemo/diatmb.f90 FFLAGS__nemo__diatmb.flags dom_oce.o in_out_manager.o iom.o oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diawri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diawri.done: diawri.o diahth.done dianam.done dom_oce.done in_out_manager.done ioipsl.done iom.done lbclnk.done lib_mpp.done limwri.done oce.done phycst.done sbc_ice.done sbc_oce.done sbcssr.done timing.done wrk_nemo.done zdf_oce.done zdfddm.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

diawri.o: $(FCM_PPSRCDIR)/nemo/diawri.f90 FFLAGS__nemo__diawri.flags diahth.o dianam.o dom_oce.o in_out_manager.o ioipsl.o iom.o lbclnk.o lib_mpp.o limwri.o oce.o phycst.o sbc_ice.o sbc_oce.o sbcssr.o timing.o wrk_nemo.o zdf_oce.o zdfddm.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__diurnal_bulk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

diurnal_bulk.done: diurnal_bulk.o dom_oce.done in_out_manager.done lib_mpp.done par_kind.done phycst.done solfrac_mod.done
	touch $(FCM_DONEDIR)/$@

diurnal_bulk.o: $(FCM_PPSRCDIR)/nemo/diurnal_bulk.f90 FFLAGS__nemo__diurnal_bulk.flags dom_oce.o in_out_manager.o lib_mpp.o par_kind.o phycst.o solfrac_mod.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__divhor.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

divhor.done: divhor.o dom_oce.done in_out_manager.done iscplhsb.done iscplini.done lbclnk.done lib_mpp.done oce.done sbc_oce.done sbcisf.done sbcrnf.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

divhor.o: $(FCM_PPSRCDIR)/nemo/divhor.f90 FFLAGS__nemo__divhor.flags dom_oce.o in_out_manager.o iscplhsb.o iscplini.o lbclnk.o lib_mpp.o oce.o sbc_oce.o sbcisf.o sbcrnf.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dom_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dom_oce.done: dom_oce.o par_oce.done
	touch $(FCM_DONEDIR)/$@

dom_oce.o: $(FCM_PPSRCDIR)/nemo/dom_oce.f90 FFLAGS__nemo__dom_oce.flags par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domain.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domain.done: domain.o c1d.done dom_oce.done domc1d.done domhgr.done dommsk.done domvvl.done domwri.done domzgr.done dyncor_c1d.done in_out_manager.done ioipsl.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done sbc_oce.done timing.done trc_oce.done usrdef_closea.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domain.o: $(FCM_PPSRCDIR)/nemo/domain.f90 FFLAGS__nemo__domain.flags c1d.o dom_oce.o domc1d.o domhgr.o dommsk.o domvvl.o domwri.o domzgr.o dyncor_c1d.o in_out_manager.o ioipsl.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o sbc_oce.o timing.o trc_oce.o usrdef_closea.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domc1d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domc1d.done: domc1d.o
	touch $(FCM_DONEDIR)/$@

domc1d.o: $(FCM_PPSRCDIR)/nemo/domc1d.f90 FFLAGS__nemo__domc1d.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domhgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domhgr.done: domhgr.o dom_oce.done in_out_manager.done iom.done lib_mpp.done par_oce.done phycst.done timing.done usrdef_hgr.done
	touch $(FCM_DONEDIR)/$@

domhgr.o: $(FCM_PPSRCDIR)/nemo/domhgr.f90 FFLAGS__nemo__domhgr.flags dom_oce.o in_out_manager.o iom.o lib_mpp.o par_oce.o phycst.o timing.o usrdef_hgr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dommsk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dommsk.done: dommsk.o bdy_oce.done dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done timing.done usrdef_fmask.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dommsk.o: $(FCM_PPSRCDIR)/nemo/dommsk.f90 FFLAGS__nemo__dommsk.flags bdy_oce.o dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o timing.o usrdef_fmask.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domngb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domngb.done: domngb.o dom_oce.done in_out_manager.done lib_mpp.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domngb.o: $(FCM_PPSRCDIR)/nemo/domngb.f90 FFLAGS__nemo__domngb.flags dom_oce.o in_out_manager.o lib_mpp.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domvvl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domvvl.done: domvvl.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done restart.done sbc_oce.done timing.done usrdef_istate.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domvvl.o: $(FCM_PPSRCDIR)/nemo/domvvl.f90 FFLAGS__nemo__domvvl.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o restart.o sbc_oce.o timing.o usrdef_istate.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domwri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domwri.done: domwri.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done phycst.done timing.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domwri.o: $(FCM_PPSRCDIR)/nemo/domwri.f90 FFLAGS__nemo__domwri.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o phycst.o timing.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domzgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domzgr.done: domzgr.o depth_e3.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done timing.done usrdef_zgr.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domzgr.o: $(FCM_PPSRCDIR)/nemo/domzgr.f90 FFLAGS__nemo__domzgr.flags depth_e3.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o timing.o usrdef_zgr.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dtatsd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dtatsd.done: dtatsd.o dom_oce.done fldread.done in_out_manager.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dtatsd.o: $(FCM_PPSRCDIR)/nemo/dtatsd.f90 FFLAGS__nemo__dtatsd.flags dom_oce.o fldread.o in_out_manager.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dtauvd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dtauvd.done: dtauvd.o dom_oce.done fldread.done in_out_manager.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dtauvd.o: $(FCM_PPSRCDIR)/nemo/dtauvd.f90 FFLAGS__nemo__dtauvd.flags dom_oce.o fldread.o in_out_manager.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynadv.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynadv.done: dynadv.o dom_oce.done dynadv_cen2.done dynadv_ubs.done dynkeg.done dynzad.done in_out_manager.done lib_mpp.done timing.done
	touch $(FCM_DONEDIR)/$@

dynadv.o: $(FCM_PPSRCDIR)/nemo/dynadv.f90 FFLAGS__nemo__dynadv.flags dom_oce.o dynadv_cen2.o dynadv_ubs.o dynkeg.o dynzad.o in_out_manager.o lib_mpp.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynadv_cen2.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynadv_cen2.done: dynadv_cen2.o dom_oce.done in_out_manager.done lib_mpp.done oce.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynadv_cen2.o: $(FCM_PPSRCDIR)/nemo/dynadv_cen2.f90 FFLAGS__nemo__dynadv_cen2.flags dom_oce.o in_out_manager.o lib_mpp.o oce.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynadv_ubs.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynadv_ubs.done: dynadv_ubs.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynadv_ubs.o: $(FCM_PPSRCDIR)/nemo/dynadv_ubs.f90 FFLAGS__nemo__dynadv_ubs.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynbfr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynbfr.done: dynbfr.o dom_oce.done in_out_manager.done oce.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done zdf_oce.done zdfbfr.done
	touch $(FCM_DONEDIR)/$@

dynbfr.o: $(FCM_PPSRCDIR)/nemo/dynbfr.f90 FFLAGS__nemo__dynbfr.flags dom_oce.o in_out_manager.o oce.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o zdf_oce.o zdfbfr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dyncor_c1d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dyncor_c1d.done: dyncor_c1d.o
	touch $(FCM_DONEDIR)/$@

dyncor_c1d.o: $(FCM_PPSRCDIR)/nemo/dyncor_c1d.f90 FFLAGS__nemo__dyncor_c1d.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dyndmp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dyndmp.done: dyndmp.o c1d.done dom_oce.done dtauvd.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done prtctl.done timing.done tradmp.done wrk_nemo.done zdf_oce.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

dyndmp.o: $(FCM_PPSRCDIR)/nemo/dyndmp.f90 FFLAGS__nemo__dyndmp.flags c1d.o dom_oce.o dtauvd.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o prtctl.o timing.o tradmp.o wrk_nemo.o zdf_oce.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynhpg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynhpg.done: dynhpg.o dom_oce.done eosbn2.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done timing.done trd_oce.done trddyn.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynhpg.o: $(FCM_PPSRCDIR)/nemo/dynhpg.f90 FFLAGS__nemo__dynhpg.flags dom_oce.o eosbn2.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o timing.o trd_oce.o trddyn.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynkeg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynkeg.done: dynkeg.o bdy_oce.done dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynkeg.o: $(FCM_PPSRCDIR)/nemo/dynkeg.f90 FFLAGS__nemo__dynkeg.flags bdy_oce.o dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynldf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynldf.done: dynldf.o dom_oce.done dynldf_iso.done dynldf_lap_blp.done in_out_manager.done lbclnk.done ldfdyn.done ldfslp.done lib_mpp.done oce.done phycst.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynldf.o: $(FCM_PPSRCDIR)/nemo/dynldf.f90 FFLAGS__nemo__dynldf.flags dom_oce.o dynldf_iso.o dynldf_lap_blp.o in_out_manager.o lbclnk.o ldfdyn.o ldfslp.o lib_mpp.o oce.o phycst.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynldf_iso.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynldf_iso.done: dynldf_iso.o dom_oce.done in_out_manager.done lbclnk.done ldfdyn.done ldfslp.done ldftra.done lib_mpp.done oce.done prtctl.done timing.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

dynldf_iso.o: $(FCM_PPSRCDIR)/nemo/dynldf_iso.f90 FFLAGS__nemo__dynldf_iso.flags dom_oce.o in_out_manager.o lbclnk.o ldfdyn.o ldfslp.o ldftra.o lib_mpp.o oce.o prtctl.o timing.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynldf_lap_blp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynldf_lap_blp.done: dynldf_lap_blp.o dom_oce.done in_out_manager.done lbclnk.done ldfdyn.done ldfslp.done oce.done timing.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

dynldf_lap_blp.o: $(FCM_PPSRCDIR)/nemo/dynldf_lap_blp.f90 FFLAGS__nemo__dynldf_lap_blp.flags dom_oce.o in_out_manager.o lbclnk.o ldfdyn.o ldfslp.o oce.o timing.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynnxt.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynnxt.done: dynnxt.o bdy_oce.done bdydta.done bdydyn.done bdyvol.done dom_oce.done domvvl.done dynadv.done dynspg_ts.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done timing.done trd_oce.done trddyn.done trdken.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynnxt.o: $(FCM_PPSRCDIR)/nemo/dynnxt.f90 FFLAGS__nemo__dynnxt.flags bdy_oce.o bdydta.o bdydyn.o bdyvol.o dom_oce.o domvvl.o dynadv.o dynspg_ts.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o timing.o trd_oce.o trddyn.o trdken.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynspg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynspg.done: dynspg.o c1d.done dom_oce.done dynspg_exp.done dynspg_ts.done in_out_manager.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done sbcapr.done sbctide.done timing.done trd_oce.done trddyn.done updtide.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynspg.o: $(FCM_PPSRCDIR)/nemo/dynspg.f90 FFLAGS__nemo__dynspg.flags c1d.o dom_oce.o dynspg_exp.o dynspg_ts.o in_out_manager.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o sbcapr.o sbctide.o timing.o trd_oce.o trddyn.o updtide.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynspg_exp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynspg_exp.done: dynspg_exp.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

dynspg_exp.o: $(FCM_PPSRCDIR)/nemo/dynspg_exp.f90 FFLAGS__nemo__dynspg_exp.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynspg_ts.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynspg_ts.done: dynspg_ts.o bdy_oce.done bdydyn2d.done bdytides.done diatmb.done dom_oce.done dynadv.done dynvor.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done restart.done sbc_oce.done sbcapr.done sbcisf.done sbctide.done sbcwave.done timing.done updtide.done wet_dry.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

dynspg_ts.o: $(FCM_PPSRCDIR)/nemo/dynspg_ts.f90 FFLAGS__nemo__dynspg_ts.flags bdy_oce.o bdydyn2d.o bdytides.o diatmb.o dom_oce.o dynadv.o dynvor.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o restart.o sbc_oce.o sbcapr.o sbcisf.o sbctide.o sbcwave.o timing.o updtide.o wet_dry.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynvor.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynvor.done: dynvor.o dom_oce.done dommsk.done dynadv.done in_out_manager.done lbclnk.done lib_mpp.done oce.done prtctl.done sbc_oce.done sbcwave.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynvor.o: $(FCM_PPSRCDIR)/nemo/dynvor.f90 FFLAGS__nemo__dynvor.flags dom_oce.o dommsk.o dynadv.o in_out_manager.o lbclnk.o lib_mpp.o oce.o prtctl.o sbc_oce.o sbcwave.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynzad.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynzad.done: dynzad.o dom_oce.done in_out_manager.done lib_mpp.done oce.done prtctl.done sbc_oce.done timing.done trd_oce.done trddyn.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dynzad.o: $(FCM_PPSRCDIR)/nemo/dynzad.f90 FFLAGS__nemo__dynzad.flags dom_oce.o in_out_manager.o lib_mpp.o oce.o prtctl.o sbc_oce.o timing.o trd_oce.o trddyn.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynzdf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynzdf.done: dynzdf.o dom_oce.done dynzdf_exp.done dynzdf_imp.done in_out_manager.done ldfdyn.done lib_mpp.done oce.done prtctl.done timing.done trd_oce.done trddyn.done wrk_nemo.done zdf_oce.done zdfgls.done zdftke.done
	touch $(FCM_DONEDIR)/$@

dynzdf.o: $(FCM_PPSRCDIR)/nemo/dynzdf.f90 FFLAGS__nemo__dynzdf.flags dom_oce.o dynzdf_exp.o dynzdf_imp.o in_out_manager.o ldfdyn.o lib_mpp.o oce.o prtctl.o timing.o trd_oce.o trddyn.o wrk_nemo.o zdf_oce.o zdfgls.o zdftke.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynzdf_exp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynzdf_exp.done: dynzdf_exp.o dom_oce.done dynadv.done in_out_manager.done lib_mpp.done oce.done phycst.done sbc_oce.done timing.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

dynzdf_exp.o: $(FCM_PPSRCDIR)/nemo/dynzdf_exp.f90 FFLAGS__nemo__dynzdf_exp.flags dom_oce.o dynadv.o in_out_manager.o lib_mpp.o oce.o phycst.o sbc_oce.o timing.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dynzdf_imp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dynzdf_imp.done: dynzdf_imp.o dom_oce.done domvvl.done dynadv.done in_out_manager.done lib_mpp.done oce.done phycst.done sbc_oce.done timing.done wrk_nemo.done zdf_oce.done zdfbfr.done
	touch $(FCM_DONEDIR)/$@

dynzdf_imp.o: $(FCM_PPSRCDIR)/nemo/dynzdf_imp.f90 FFLAGS__nemo__dynzdf_imp.flags dom_oce.o domvvl.o dynadv.o in_out_manager.o lib_mpp.o oce.o phycst.o sbc_oce.o timing.o wrk_nemo.o zdf_oce.o zdfbfr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__eosbn2.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

eosbn2.done: eosbn2.o dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done phycst.done prtctl.done stopar.done stopts.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

eosbn2.o: $(FCM_PPSRCDIR)/nemo/eosbn2.f90 FFLAGS__nemo__eosbn2.flags dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o phycst.o prtctl.o stopar.o stopts.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

find_obs_proc.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/find_obs_proc.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

find_obs_proc.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/find_obs_proc.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__fldread.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

fldread.done: fldread.o bdy_oce.done dom_oce.done geo2ocean.done in_out_manager.done ioipsl.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

fldread.o: $(FCM_PPSRCDIR)/nemo/fldread.f90 FFLAGS__nemo__fldread.flags bdy_oce.o dom_oce.o geo2ocean.o in_out_manager.o ioipsl.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__flo4rk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

flo4rk.done: flo4rk.o
	touch $(FCM_DONEDIR)/$@

flo4rk.o: $(FCM_PPSRCDIR)/nemo/flo4rk.f90 FFLAGS__nemo__flo4rk.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__flo_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

flo_oce.done: flo_oce.o
	touch $(FCM_DONEDIR)/$@

flo_oce.o: $(FCM_PPSRCDIR)/nemo/flo_oce.f90 FFLAGS__nemo__flo_oce.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__floats.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

floats.done: floats.o
	touch $(FCM_DONEDIR)/$@

floats.o: $(FCM_PPSRCDIR)/nemo/floats.f90 FFLAGS__nemo__floats.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__floblk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

floblk.done: floblk.o
	touch $(FCM_DONEDIR)/$@

floblk.o: $(FCM_PPSRCDIR)/nemo/floblk.f90 FFLAGS__nemo__floblk.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__flodom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

flodom.done: flodom.o
	touch $(FCM_DONEDIR)/$@

flodom.o: $(FCM_PPSRCDIR)/nemo/flodom.f90 FFLAGS__nemo__flodom.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__florst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

florst.done: florst.o
	touch $(FCM_DONEDIR)/$@

florst.o: $(FCM_PPSRCDIR)/nemo/florst.f90 FFLAGS__nemo__florst.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__flowri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

flowri.done: flowri.o
	touch $(FCM_DONEDIR)/$@

flowri.o: $(FCM_PPSRCDIR)/nemo/flowri.f90 FFLAGS__nemo__flowri.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__geo2ocean.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

geo2ocean.done: geo2ocean.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done phycst.done
	touch $(FCM_DONEDIR)/$@

geo2ocean.o: $(FCM_PPSRCDIR)/nemo/geo2ocean.f90 FFLAGS__nemo__geo2ocean.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o phycst.o
	fcm_internal compile:F nemo $< $@

greg2jul.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/greg2jul.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

greg2jul.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/greg2jul.h90
	touch $(FCM_DONEDIR)/$@

grt_cir_dis.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/grt_cir_dis.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

grt_cir_dis.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/grt_cir_dis.h90
	touch $(FCM_DONEDIR)/$@

grt_cir_dis_saa.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/grt_cir_dis_saa.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

grt_cir_dis_saa.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/grt_cir_dis_saa.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__icb_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icb_oce.done: icb_oce.o lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

icb_oce.o: $(FCM_PPSRCDIR)/nemo/icb_oce.f90 FFLAGS__nemo__icb_oce.flags lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbclv.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbclv.done: icbclv.o dom_oce.done icb_oce.done icbdia.done icbutl.done lbclnk.done lib_mpp.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

icbclv.o: $(FCM_PPSRCDIR)/nemo/icbclv.f90 FFLAGS__nemo__icbclv.flags dom_oce.o icb_oce.o icbdia.o icbutl.o lbclnk.o lib_mpp.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbdia.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbdia.done: icbdia.o dom_oce.done icb_oce.done icbutl.done in_out_manager.done iom.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

icbdia.o: $(FCM_PPSRCDIR)/nemo/icbdia.f90 FFLAGS__nemo__icbdia.flags dom_oce.o icb_oce.o icbutl.o in_out_manager.o iom.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbdyn.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbdyn.done: icbdyn.o dom_oce.done icb_oce.done icbdia.done icbutl.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

icbdyn.o: $(FCM_PPSRCDIR)/nemo/icbdyn.f90 FFLAGS__nemo__icbdyn.flags dom_oce.o icb_oce.o icbdia.o icbutl.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbini.done: icbini.o dom_oce.done fldread.done icb_oce.done icbdia.done icbrst.done icbtrj.done icbutl.done in_out_manager.done iom.done lbclnk.done lib_mpp.done sbc_ice.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

icbini.o: $(FCM_PPSRCDIR)/nemo/icbini.f90 FFLAGS__nemo__icbini.flags dom_oce.o fldread.o icb_oce.o icbdia.o icbrst.o icbtrj.o icbutl.o in_out_manager.o iom.o lbclnk.o lib_mpp.o sbc_ice.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icblbc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icblbc.done: icblbc.o dom_oce.done icb_oce.done icbutl.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

icblbc.o: $(FCM_PPSRCDIR)/nemo/icblbc.f90 FFLAGS__nemo__icblbc.flags dom_oce.o icb_oce.o icbutl.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbrst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbrst.done: icbrst.o dom_oce.done icb_oce.done icbutl.done in_out_manager.done iom.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

icbrst.o: $(FCM_PPSRCDIR)/nemo/icbrst.f90 FFLAGS__nemo__icbrst.flags dom_oce.o icb_oce.o icbutl.o in_out_manager.o iom.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbstp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbstp.done: icbstp.o dom_oce.done fldread.done icb_oce.done icbclv.done icbdia.done icbdyn.done icbini.done icblbc.done icbrst.done icbthm.done icbtrj.done icbutl.done in_out_manager.done iom.done lib_mpp.done par_oce.done phycst.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

icbstp.o: $(FCM_PPSRCDIR)/nemo/icbstp.f90 FFLAGS__nemo__icbstp.flags dom_oce.o fldread.o icb_oce.o icbclv.o icbdia.o icbdyn.o icbini.o icblbc.o icbrst.o icbthm.o icbtrj.o icbutl.o in_out_manager.o iom.o lib_mpp.o par_oce.o phycst.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbthm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbthm.done: icbthm.o dom_oce.done icb_oce.done icbdia.done icbutl.done in_out_manager.done lib_mpp.done par_oce.done phycst.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

icbthm.o: $(FCM_PPSRCDIR)/nemo/icbthm.f90 FFLAGS__nemo__icbthm.flags dom_oce.o icb_oce.o icbdia.o icbutl.o in_out_manager.o lib_mpp.o par_oce.o phycst.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbtrj.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbtrj.done: icbtrj.o dom_oce.done icb_oce.done icbutl.done in_out_manager.done lib_mpp.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

icbtrj.o: $(FCM_PPSRCDIR)/nemo/icbtrj.f90 FFLAGS__nemo__icbtrj.flags dom_oce.o icb_oce.o icbutl.o in_out_manager.o lib_mpp.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__icbutl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

icbutl.done: icbutl.o dom_oce.done icb_oce.done ice.done in_out_manager.done lbclnk.done lib_mpp.done par_oce.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

icbutl.o: $(FCM_PPSRCDIR)/nemo/icbutl.f90 FFLAGS__nemo__icbutl.flags dom_oce.o icb_oce.o ice.o in_out_manager.o lbclnk.o lib_mpp.o par_oce.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ice.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ice.done: ice.o in_out_manager.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

ice.o: $(FCM_PPSRCDIR)/nemo/ice.f90 FFLAGS__nemo__ice.flags in_out_manager.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__in_out_manager.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

in_out_manager.done: in_out_manager.o nc4interface.done par_oce.done
	touch $(FCM_DONEDIR)/$@

in_out_manager.o: $(FCM_PPSRCDIR)/nemo/in_out_manager.f90 FFLAGS__nemo__in_out_manager.flags nc4interface.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom.done: iom.o c1d.done crs.done dianam.done dom_oce.done domngb.done flo_oce.done icb_oce.done ice.done in_out_manager.done ioipsl.done iom_def.done iom_nf90.done lbclnk.done lib_mpp.done phycst.done sbc_oce.done trc_oce.done
	touch $(FCM_DONEDIR)/$@

iom.o: $(FCM_PPSRCDIR)/nemo/iom.f90 FFLAGS__nemo__iom.flags c1d.o crs.o dianam.o dom_oce.o domngb.o flo_oce.o icb_oce.o ice.o in_out_manager.o ioipsl.o iom_def.o iom_nf90.o lbclnk.o lib_mpp.o phycst.o sbc_oce.o trc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom_def.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom_def.done: iom_def.o par_kind.done
	touch $(FCM_DONEDIR)/$@

iom_def.o: $(FCM_PPSRCDIR)/nemo/iom_def.f90 FFLAGS__nemo__iom_def.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom_nf90.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom_nf90.done: iom_nf90.o dom_oce.done in_out_manager.done iom_def.done lbclnk.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

iom_nf90.o: $(FCM_PPSRCDIR)/nemo/iom_nf90.f90 FFLAGS__nemo__iom_nf90.flags dom_oce.o in_out_manager.o iom_def.o lbclnk.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplhsb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplhsb.done: iscplhsb.o dom_oce.done domngb.done domwri.done in_out_manager.done iscplini.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

iscplhsb.o: $(FCM_PPSRCDIR)/nemo/iscplhsb.f90 FFLAGS__nemo__iscplhsb.flags dom_oce.o domngb.o domwri.o in_out_manager.o iscplini.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplini.done: iscplini.o dom_oce.done in_out_manager.done lib_fortran.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

iscplini.o: $(FCM_PPSRCDIR)/nemo/iscplini.f90 FFLAGS__nemo__iscplini.flags dom_oce.o in_out_manager.o lib_fortran.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplrst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplrst.done: iscplrst.o dom_oce.done domvvl.done domwri.done in_out_manager.done iom.done iscplhsb.done iscplini.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

iscplrst.o: $(FCM_PPSRCDIR)/nemo/iscplrst.f90 FFLAGS__nemo__iscplrst.flags dom_oce.o domvvl.o domwri.o in_out_manager.o iom.o iscplhsb.o iscplini.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__istate.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

istate.done: istate.o daymod.done divhor.done dom_oce.done domvvl.done dtatsd.done dtauvd.done in_out_manager.done iom.done iscplrst.done lib_mpp.done oce.done restart.done timing.done usrdef_istate.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

istate.o: $(FCM_PPSRCDIR)/nemo/istate.f90 FFLAGS__nemo__istate.flags daymod.o divhor.o dom_oce.o domvvl.o dtatsd.o dtauvd.o in_out_manager.o iom.o iscplrst.o lib_mpp.o oce.o restart.o timing.o usrdef_istate.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

jul2greg.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/jul2greg.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

jul2greg.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/jul2greg.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__julian.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

julian.done: julian.o lib_mpp.done par_kind.done
	touch $(FCM_DONEDIR)/$@

julian.o: $(FCM_PPSRCDIR)/nemo/julian.f90 FFLAGS__nemo__julian.flags lib_mpp.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lbclnk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lbclnk.done: lbclnk.o lib_mpp.done
	touch $(FCM_DONEDIR)/$@

lbclnk.o: $(FCM_PPSRCDIR)/nemo/lbclnk.f90 FFLAGS__nemo__lbclnk.flags lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lbcnfd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lbcnfd.done: lbcnfd.o dom_oce.done in_out_manager.done
	touch $(FCM_DONEDIR)/$@

lbcnfd.o: $(FCM_PPSRCDIR)/nemo/lbcnfd.f90 FFLAGS__nemo__lbcnfd.flags dom_oce.o in_out_manager.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ldfc1d_c2d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ldfc1d_c2d.done: ldfc1d_c2d.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

ldfc1d_c2d.o: $(FCM_PPSRCDIR)/nemo/ldfc1d_c2d.f90 FFLAGS__nemo__ldfc1d_c2d.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ldfdyn.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ldfdyn.done: ldfdyn.o dom_oce.done in_out_manager.done iom.done lbclnk.done ldfc1d_c2d.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

ldfdyn.o: $(FCM_PPSRCDIR)/nemo/ldfdyn.f90 FFLAGS__nemo__ldfdyn.flags dom_oce.o in_out_manager.o iom.o lbclnk.o ldfc1d_c2d.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ldfslp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ldfslp.done: ldfslp.o dom_oce.done eosbn2.done in_out_manager.done lbclnk.done ldfdyn.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done timing.done wrk_nemo.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

ldfslp.o: $(FCM_PPSRCDIR)/nemo/ldfslp.f90 FFLAGS__nemo__ldfslp.flags dom_oce.o eosbn2.o in_out_manager.o lbclnk.o ldfdyn.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o timing.o wrk_nemo.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ldftra.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ldftra.done: ldftra.o diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done ldfc1d_c2d.done ldfslp.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

ldftra.o: $(FCM_PPSRCDIR)/nemo/ldftra.f90 FFLAGS__nemo__ldftra.flags diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o ldfc1d_c2d.o ldfslp.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_cray.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_cray.done: lib_cray.o
	touch $(FCM_DONEDIR)/$@

lib_cray.o: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/lib_cray.f90 FFLAGS__nemo__lib_cray.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_fortran.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_fortran.done: lib_fortran.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

lib_fortran.o: $(FCM_PPSRCDIR)/nemo/lib_fortran.f90 FFLAGS__nemo__lib_fortran.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_mpp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_mpp.done: lib_mpp.o dom_oce.done in_out_manager.done lbcnfd.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

lib_mpp.o: $(FCM_PPSRCDIR)/nemo/lib_mpp.f90 FFLAGS__nemo__lib_mpp.flags dom_oce.o in_out_manager.o lbcnfd.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limadv_prather.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limadv_prather.done: limadv_prather.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done prtctl.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limadv_prather.o: $(FCM_PPSRCDIR)/nemo/limadv_prather.f90 FFLAGS__nemo__limadv_prather.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o prtctl.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limadv_umx.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limadv_umx.done: limadv_umx.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done phycst.done sbc_oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limadv_umx.o: $(FCM_PPSRCDIR)/nemo/limadv_umx.f90 FFLAGS__nemo__limadv_umx.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o phycst.o sbc_oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limcons.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limcons.done: limcons.o dom_oce.done ice.done in_out_manager.done lib_fortran.done lib_mpp.done phycst.done sbc_ice.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

limcons.o: $(FCM_PPSRCDIR)/nemo/limcons.f90 FFLAGS__nemo__limcons.flags dom_oce.o ice.o in_out_manager.o lib_fortran.o lib_mpp.o phycst.o sbc_ice.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limctl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limctl.done: limctl.o dom_oce.done ice.done in_out_manager.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done thd_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limctl.o: $(FCM_PPSRCDIR)/nemo/limctl.f90 FFLAGS__nemo__limctl.flags dom_oce.o ice.o in_out_manager.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o thd_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limdiahsb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limdiahsb.done: limdiahsb.o daymod.done dom_oce.done ice.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done limrst.done phycst.done sbc_ice.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

limdiahsb.o: $(FCM_PPSRCDIR)/nemo/limdiahsb.f90 FFLAGS__nemo__limdiahsb.flags daymod.o dom_oce.o ice.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o limrst.o phycst.o sbc_ice.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limdyn.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limdyn.done: limdyn.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limctl.done limrhg.done limvar.done phycst.done sbc_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limdyn.o: $(FCM_PPSRCDIR)/nemo/limdyn.f90 FFLAGS__nemo__limdyn.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limctl.o limrhg.o limvar.o phycst.o sbc_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limhdf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limhdf.done: limhdf.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done prtctl.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limhdf.o: $(FCM_PPSRCDIR)/nemo/limhdf.f90 FFLAGS__nemo__limhdf.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o prtctl.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limistate.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limistate.done: limistate.o dom_oce.done eosbn2.done fldread.done ice.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done limvar.done oce.done par_oce.done phycst.done sbc_ice.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limistate.o: $(FCM_PPSRCDIR)/nemo/limistate.f90 FFLAGS__nemo__limistate.flags dom_oce.o eosbn2.o fldread.o ice.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o limvar.o oce.o par_oce.o phycst.o sbc_ice.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limitd_me.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limitd_me.done: limitd_me.o dom_oce.done ice.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limctl.done limvar.done par_oce.done phycst.done sbc_oce.done thd_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limitd_me.o: $(FCM_PPSRCDIR)/nemo/limitd_me.f90 FFLAGS__nemo__limitd_me.flags dom_oce.o ice.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limctl.o limvar.o par_oce.o phycst.o sbc_oce.o thd_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limitd_th.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limitd_th.done: limitd_th.o dom_oce.done ice.done in_out_manager.done lib_fortran.done lib_mpp.done limcons.done limvar.done par_oce.done phycst.done prtctl.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limitd_th.o: $(FCM_PPSRCDIR)/nemo/limitd_th.f90 FFLAGS__nemo__limitd_th.flags dom_oce.o ice.o in_out_manager.o lib_fortran.o lib_mpp.o limcons.o limvar.o par_oce.o phycst.o prtctl.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limrhg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limrhg.done: limrhg.o bdy_oce.done bdyice_lim.done dom_oce.done ice.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limitd_me.done oce.done par_oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limrhg.o: $(FCM_PPSRCDIR)/nemo/limrhg.f90 FFLAGS__nemo__limrhg.flags bdy_oce.o bdyice_lim.o dom_oce.o ice.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limitd_me.o oce.o par_oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limrst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limrst.done: limrst.o dom_oce.done ice.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done limctl.done oce.done sbc_ice.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limrst.o: $(FCM_PPSRCDIR)/nemo/limrst.f90 FFLAGS__nemo__limrst.flags dom_oce.o ice.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o limctl.o oce.o sbc_ice.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limsbc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limsbc.done: limsbc.o albedo.done bdy_oce.done dom_oce.done domvvl.done ice.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limctl.done oce.done par_oce.done phycst.done sbc_ice.done sbc_oce.done sbccpl.done traqsr.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limsbc.o: $(FCM_PPSRCDIR)/nemo/limsbc.f90 FFLAGS__nemo__limsbc.flags albedo.o bdy_oce.o dom_oce.o domvvl.o ice.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limctl.o oce.o par_oce.o phycst.o sbc_ice.o sbc_oce.o sbccpl.o traqsr.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limtab.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limtab.done: limtab.o par_kind.done
	touch $(FCM_DONEDIR)/$@

limtab.o: $(FCM_PPSRCDIR)/nemo/limtab.f90 FFLAGS__nemo__limtab.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd.done: limthd.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limctl.done limitd_th.done limtab.done limthd_da.done limthd_dh.done limthd_dif.done limthd_ent.done limthd_lac.done limthd_sal.done limvar.done phycst.done sbc_ice.done sbc_oce.done thd_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd.o: $(FCM_PPSRCDIR)/nemo/limthd.f90 FFLAGS__nemo__limthd.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limctl.o limitd_th.o limtab.o limthd_da.o limthd_dh.o limthd_dif.o limthd_ent.o limthd_lac.o limthd_sal.o limvar.o phycst.o sbc_ice.o sbc_oce.o thd_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_da.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_da.done: limthd_da.o ice.done lib_fortran.done lib_mpp.done par_oce.done phycst.done sbc_ice.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_da.o: $(FCM_PPSRCDIR)/nemo/limthd_da.f90 FFLAGS__nemo__limthd_da.flags ice.o lib_fortran.o lib_mpp.o par_oce.o phycst.o sbc_ice.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_dh.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_dh.done: limthd_dh.o ice.done in_out_manager.done lib_fortran.done lib_mpp.done par_oce.done phycst.done sbc_oce.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_dh.o: $(FCM_PPSRCDIR)/nemo/limthd_dh.f90 FFLAGS__nemo__limthd_dh.flags ice.o in_out_manager.o lib_fortran.o lib_mpp.o par_oce.o phycst.o sbc_oce.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_dif.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_dif.done: limthd_dif.o ice.done in_out_manager.done lib_fortran.done lib_mpp.done par_oce.done phycst.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_dif.o: $(FCM_PPSRCDIR)/nemo/limthd_dif.f90 FFLAGS__nemo__limthd_dif.flags ice.o in_out_manager.o lib_fortran.o lib_mpp.o par_oce.o phycst.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_ent.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_ent.done: limthd_ent.o dom_oce.done domain.done ice.done in_out_manager.done lib_fortran.done lib_mpp.done limvar.done par_oce.done phycst.done sbc_oce.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_ent.o: $(FCM_PPSRCDIR)/nemo/limthd_ent.f90 FFLAGS__nemo__limthd_ent.flags dom_oce.o domain.o ice.o in_out_manager.o lib_fortran.o lib_mpp.o limvar.o par_oce.o phycst.o sbc_oce.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_lac.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_lac.done: limthd_lac.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limtab.done limthd_ent.done limvar.done par_oce.done phycst.done sbc_ice.done sbc_oce.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_lac.o: $(FCM_PPSRCDIR)/nemo/limthd_lac.f90 FFLAGS__nemo__limthd_lac.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limtab.o limthd_ent.o limvar.o par_oce.o phycst.o sbc_ice.o sbc_oce.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limthd_sal.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limthd_sal.done: limthd_sal.o ice.done in_out_manager.done lib_fortran.done lib_mpp.done limvar.done par_oce.done phycst.done sbc_oce.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limthd_sal.o: $(FCM_PPSRCDIR)/nemo/limthd_sal.f90 FFLAGS__nemo__limthd_sal.flags ice.o in_out_manager.o lib_fortran.o lib_mpp.o limvar.o par_oce.o phycst.o sbc_oce.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limtrp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limtrp.done: limtrp.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limadv_prather.done limadv_umx.done limcons.done limctl.done limhdf.done limvar.done phycst.done prtctl.done sbc_oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limtrp.o: $(FCM_PPSRCDIR)/nemo/limtrp.f90 FFLAGS__nemo__limtrp.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limadv_prather.o limadv_umx.o limcons.o limctl.o limhdf.o limvar.o phycst.o prtctl.o sbc_oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limupdate1.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limupdate1.done: limupdate1.o dom_oce.done ice.done in_out_manager.done lib_fortran.done lib_mpp.done limcons.done limctl.done limitd_th.done limvar.done phycst.done sbc_ice.done sbc_oce.done thd_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limupdate1.o: $(FCM_PPSRCDIR)/nemo/limupdate1.f90 FFLAGS__nemo__limupdate1.flags dom_oce.o ice.o in_out_manager.o lib_fortran.o lib_mpp.o limcons.o limctl.o limitd_th.o limvar.o phycst.o sbc_ice.o sbc_oce.o thd_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limupdate2.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limupdate2.done: limupdate2.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limcons.done limctl.done limitd_th.done limvar.done phycst.done sbc_ice.done sbc_oce.done thd_ice.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limupdate2.o: $(FCM_PPSRCDIR)/nemo/limupdate2.f90 FFLAGS__nemo__limupdate2.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limcons.o limctl.o limitd_th.o limvar.o phycst.o sbc_ice.o sbc_oce.o thd_ice.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limvar.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limvar.done: limvar.o ice.done in_out_manager.done lib_fortran.done lib_mpp.done par_oce.done phycst.done sbc_ice.done sbc_oce.done thd_ice.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limvar.o: $(FCM_PPSRCDIR)/nemo/limvar.f90 FFLAGS__nemo__limvar.flags ice.o in_out_manager.o lib_fortran.o lib_mpp.o par_oce.o phycst.o sbc_ice.o sbc_oce.o thd_ice.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__limwri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

limwri.done: limwri.o dianam.done dom_oce.done ice.done in_out_manager.done ioipsl.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limvar.done phycst.done sbc_ice.done sbc_oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

limwri.o: $(FCM_PPSRCDIR)/nemo/limwri.f90 FFLAGS__nemo__limwri.flags dianam.o dom_oce.o ice.o in_out_manager.o ioipsl.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limvar.o phycst.o sbc_ice.o sbc_oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

linquad.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/linquad.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

linquad.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/linquad.h90
	touch $(FCM_DONEDIR)/$@

maxdist.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/maxdist.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

maxdist.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/maxdist.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__mpp_map.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

mpp_map.done: mpp_map.o dom_oce.done in_out_manager.done lib_mpp.done par_kind.done par_oce.done
	touch $(FCM_DONEDIR)/$@

mpp_map.o: $(FCM_PPSRCDIR)/nemo/mpp_map.f90 FFLAGS__nemo__mpp_map.flags dom_oce.o in_out_manager.o lib_mpp.o par_kind.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__mppini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

mppini.done: mppini.o bdy_oce.done dom_oce.done in_out_manager.done ioipsl.done iom.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

mppini.o: $(FCM_PPSRCDIR)/nemo/mppini.f90 FFLAGS__nemo__mppini.flags bdy_oce.o dom_oce.o in_out_manager.o ioipsl.o iom.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

mppini_2.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/mppini_2.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

mppini_2.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/mppini_2.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__nemo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo__nemo.flags: LDFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo__nemo.flags: LD__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

nemo.exe: nemo.o LD__nemo__nemo.flags LDFLAGS__nemo__nemo.flags $(OBJECTS) nemogcm.done
	fcm_internal load:F nemo $< $@

nemo.o: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/nemo.f90 FFLAGS__nemo__nemo.flags nemogcm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__nemogcm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

nemogcm.done: nemogcm.o bdy_oce.done bdydta.done bdyini.done cpl_oasis3.done daymod.done diawri.done dom_oce.done domain.done icbini.done icbstp.done lbcnfd.done lib_fortran.done lib_mpp.done mppini.done oce.done phycst.done sbc_oce.done sbcssm.done step.done step_oce.done usrdef_nam.done
	touch $(FCM_DONEDIR)/$@

nemogcm.o: $(FCM_PPSRCDIR)/nemo/nemogcm.f90 FFLAGS__nemo__nemogcm.flags bdy_oce.o bdydta.o bdyini.o cpl_oasis3.o daymod.o diawri.o dom_oce.o domain.o icbini.o icbstp.o lbcnfd.o lib_fortran.o lib_mpp.o mppini.o oce.o phycst.o sbc_oce.o sbcssm.o step.o step_oce.o usrdef_nam.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_const.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_const.done: obs_const.o par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_const.o: $(FCM_PPSRCDIR)/nemo/obs_const.f90 FFLAGS__nemo__obs_const.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_conv.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_conv.done: obs_conv.o par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_conv.o: $(FCM_PPSRCDIR)/nemo/obs_conv.f90 FFLAGS__nemo__obs_conv.flags par_kind.o
	fcm_internal compile:F nemo $< $@

obs_conv_functions.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_conv_functions.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

obs_conv_functions.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_conv_functions.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__obs_fbm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_fbm.done: obs_fbm.o obs_utils.done
	touch $(FCM_DONEDIR)/$@

obs_fbm.o: $(FCM_PPSRCDIR)/nemo/obs_fbm.f90 FFLAGS__nemo__obs_fbm.flags obs_utils.o
	fcm_internal compile:F nemo $< $@

obs_grd_bruteforce.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_grd_bruteforce.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

obs_grd_bruteforce.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_grd_bruteforce.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__obs_grid.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_grid.done: obs_grid.o dom_oce.done in_out_manager.done lib_mpp.done obs_const.done obs_mpp.done obs_utils.done par_kind.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

obs_grid.o: $(FCM_PPSRCDIR)/nemo/obs_grid.f90 FFLAGS__nemo__obs_grid.flags dom_oce.o in_out_manager.o lib_mpp.o obs_const.o obs_mpp.o obs_utils.o par_kind.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_inter_h2d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_inter_h2d.done: obs_inter_h2d.o in_out_manager.done lib_mpp.done obs_const.done obs_utils.done par_kind.done phycst.done
	touch $(FCM_DONEDIR)/$@

obs_inter_h2d.o: $(FCM_PPSRCDIR)/nemo/obs_inter_h2d.f90 FFLAGS__nemo__obs_inter_h2d.flags in_out_manager.o lib_mpp.o obs_const.o obs_utils.o par_kind.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_inter_sup.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_inter_sup.done: obs_inter_sup.o dom_oce.done in_out_manager.done lib_mpp.done mpp_map.done obs_grid.done obs_mpp.done par_kind.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

obs_inter_sup.o: $(FCM_PPSRCDIR)/nemo/obs_inter_sup.f90 FFLAGS__nemo__obs_inter_sup.flags dom_oce.o in_out_manager.o lib_mpp.o mpp_map.o obs_grid.o obs_mpp.o par_kind.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_inter_z1d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_inter_z1d.done: obs_inter_z1d.o par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_inter_z1d.o: $(FCM_PPSRCDIR)/nemo/obs_inter_z1d.f90 FFLAGS__nemo__obs_inter_z1d.flags par_kind.o
	fcm_internal compile:F nemo $< $@

obs_level_search.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_level_search.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

obs_level_search.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obs_level_search.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__obs_mpp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_mpp.done: obs_mpp.o dom_oce.done in_out_manager.done lib_mpp.done mpp_map.done
	touch $(FCM_DONEDIR)/$@

obs_mpp.o: $(FCM_PPSRCDIR)/nemo/obs_mpp.f90 FFLAGS__nemo__obs_mpp.flags dom_oce.o in_out_manager.o lib_mpp.o mpp_map.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_oper.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_oper.done: obs_oper.o dom_oce.done in_out_manager.done lib_mpp.done obs_const.done obs_grid.done obs_inter_h2d.done obs_inter_sup.done obs_inter_z1d.done obs_profiles_def.done obs_surf_def.done par_kind.done sbcdcy.done
	touch $(FCM_DONEDIR)/$@

obs_oper.o: $(FCM_PPSRCDIR)/nemo/obs_oper.f90 FFLAGS__nemo__obs_oper.flags dom_oce.o in_out_manager.o lib_mpp.o obs_const.o obs_grid.o obs_inter_h2d.o obs_inter_sup.o obs_inter_z1d.o obs_profiles_def.o obs_surf_def.o par_kind.o sbcdcy.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_prep.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_prep.done: obs_prep.o dom_oce.done in_out_manager.done lib_mpp.done obs_inter_sup.done obs_mpp.done obs_oper.done obs_profiles_def.done obs_surf_def.done par_kind.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

obs_prep.o: $(FCM_PPSRCDIR)/nemo/obs_prep.f90 FFLAGS__nemo__obs_prep.flags dom_oce.o in_out_manager.o lib_mpp.o obs_inter_sup.o obs_mpp.o obs_oper.o obs_profiles_def.o obs_surf_def.o par_kind.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_profiles.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_profiles.done: obs_profiles.o obs_profiles_def.done
	touch $(FCM_DONEDIR)/$@

obs_profiles.o: $(FCM_PPSRCDIR)/nemo/obs_profiles.f90 FFLAGS__nemo__obs_profiles.flags obs_profiles_def.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_profiles_def.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_profiles_def.done: obs_profiles_def.o in_out_manager.done lib_mpp.done obs_fbm.done obs_mpp.done par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_profiles_def.o: $(FCM_PPSRCDIR)/nemo/obs_profiles_def.f90 FFLAGS__nemo__obs_profiles_def.flags in_out_manager.o lib_mpp.o obs_fbm.o obs_mpp.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_read_altbias.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_read_altbias.done: obs_read_altbias.o dom_oce.done in_out_manager.done iom.done obs_inter_h2d.done obs_inter_sup.done obs_surf_def.done obs_utils.done oce.done par_kind.done par_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

obs_read_altbias.o: $(FCM_PPSRCDIR)/nemo/obs_read_altbias.f90 FFLAGS__nemo__obs_read_altbias.flags dom_oce.o in_out_manager.o iom.o obs_inter_h2d.o obs_inter_sup.o obs_surf_def.o obs_utils.o oce.o par_kind.o par_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_read_prof.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_read_prof.done: obs_read_prof.o dom_oce.done in_out_manager.done julian.done lib_mpp.done obs_conv.done obs_fbm.done obs_grid.done obs_mpp.done obs_oper.done obs_prep.done obs_profiles_def.done obs_sort.done obs_types.done obs_utils.done par_kind.done par_oce.done
	touch $(FCM_DONEDIR)/$@

obs_read_prof.o: $(FCM_PPSRCDIR)/nemo/obs_read_prof.f90 FFLAGS__nemo__obs_read_prof.flags dom_oce.o in_out_manager.o julian.o lib_mpp.o obs_conv.o obs_fbm.o obs_grid.o obs_mpp.o obs_oper.o obs_prep.o obs_profiles_def.o obs_sort.o obs_types.o obs_utils.o par_kind.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_read_surf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_read_surf.done: obs_read_surf.o dom_oce.done in_out_manager.done julian.done obs_fbm.done obs_grid.done obs_mpp.done obs_sort.done obs_surf_def.done obs_types.done obs_utils.done par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_read_surf.o: $(FCM_PPSRCDIR)/nemo/obs_read_surf.f90 FFLAGS__nemo__obs_read_surf.flags dom_oce.o in_out_manager.o julian.o obs_fbm.o obs_grid.o obs_mpp.o obs_sort.o obs_surf_def.o obs_types.o obs_utils.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_readmdt.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_readmdt.done: obs_readmdt.o dom_oce.done in_out_manager.done iom.done iom_nf90.done lib_mpp.done obs_const.done obs_inter_h2d.done obs_inter_sup.done obs_surf_def.done obs_utils.done oce.done par_kind.done par_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

obs_readmdt.o: $(FCM_PPSRCDIR)/nemo/obs_readmdt.f90 FFLAGS__nemo__obs_readmdt.flags dom_oce.o in_out_manager.o iom.o iom_nf90.o lib_mpp.o obs_const.o obs_inter_h2d.o obs_inter_sup.o obs_surf_def.o obs_utils.o oce.o par_kind.o par_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_rot_vel.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_rot_vel.done: obs_rot_vel.o dom_oce.done geo2ocean.done in_out_manager.done obs_fbm.done obs_grid.done obs_inter_h2d.done obs_inter_sup.done obs_profiles_def.done obs_utils.done par_kind.done par_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

obs_rot_vel.o: $(FCM_PPSRCDIR)/nemo/obs_rot_vel.f90 FFLAGS__nemo__obs_rot_vel.flags dom_oce.o geo2ocean.o in_out_manager.o obs_fbm.o obs_grid.o obs_inter_h2d.o obs_inter_sup.o obs_profiles_def.o obs_utils.o par_kind.o par_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_sort.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_sort.done: obs_sort.o par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_sort.o: $(FCM_PPSRCDIR)/nemo/obs_sort.f90 FFLAGS__nemo__obs_sort.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_sstbias.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_sstbias.done: obs_sstbias.o dom_oce.done in_out_manager.done iom.done obs_inter_h2d.done obs_inter_sup.done obs_surf_def.done obs_utils.done oce.done par_kind.done par_oce.done
	touch $(FCM_DONEDIR)/$@

obs_sstbias.o: $(FCM_PPSRCDIR)/nemo/obs_sstbias.f90 FFLAGS__nemo__obs_sstbias.flags dom_oce.o in_out_manager.o iom.o obs_inter_h2d.o obs_inter_sup.o obs_surf_def.o obs_utils.o oce.o par_kind.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_surf_def.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_surf_def.done: obs_surf_def.o obs_mpp.done par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_surf_def.o: $(FCM_PPSRCDIR)/nemo/obs_surf_def.f90 FFLAGS__nemo__obs_surf_def.flags obs_mpp.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_types.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_types.done: obs_types.o
	touch $(FCM_DONEDIR)/$@

obs_types.o: $(FCM_PPSRCDIR)/nemo/obs_types.f90 FFLAGS__nemo__obs_types.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_utils.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_utils.done: obs_utils.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

obs_utils.o: $(FCM_PPSRCDIR)/nemo/obs_utils.f90 FFLAGS__nemo__obs_utils.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__obs_write.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

obs_write.done: obs_write.o dom_oce.done in_out_manager.done julian.done lib_mpp.done obs_const.done obs_conv.done obs_fbm.done obs_grid.done obs_mpp.done obs_profiles_def.done obs_surf_def.done obs_types.done obs_utils.done par_kind.done
	touch $(FCM_DONEDIR)/$@

obs_write.o: $(FCM_PPSRCDIR)/nemo/obs_write.f90 FFLAGS__nemo__obs_write.flags dom_oce.o in_out_manager.o julian.o lib_mpp.o obs_const.o obs_conv.o obs_fbm.o obs_grid.o obs_mpp.o obs_profiles_def.o obs_surf_def.o obs_types.o obs_utils.o par_kind.o
	fcm_internal compile:F nemo $< $@

obsinter_h2d.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obsinter_h2d.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

obsinter_h2d.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obsinter_h2d.h90
	touch $(FCM_DONEDIR)/$@

obsinter_z1d.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obsinter_z1d.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

obsinter_z1d.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/obsinter_z1d.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

oce.done: oce.o lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

oce.o: $(FCM_PPSRCDIR)/nemo/oce.f90 FFLAGS__nemo__oce.flags lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__par_kind.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

par_kind.done: par_kind.o
	touch $(FCM_DONEDIR)/$@

par_kind.o: $(FCM_PPSRCDIR)/nemo/par_kind.f90 FFLAGS__nemo__par_kind.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__par_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

par_oce.done: par_oce.o par_kind.done
	touch $(FCM_DONEDIR)/$@

par_oce.o: $(FCM_PPSRCDIR)/nemo/par_oce.f90 FFLAGS__nemo__par_oce.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__phycst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

phycst.done: phycst.o in_out_manager.done par_oce.done
	touch $(FCM_DONEDIR)/$@

phycst.o: $(FCM_PPSRCDIR)/nemo/phycst.f90 FFLAGS__nemo__phycst.flags in_out_manager.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__prtctl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

prtctl.done: prtctl.o dom_oce.done in_out_manager.done lib_mpp.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

prtctl.o: $(FCM_PPSRCDIR)/nemo/prtctl.f90 FFLAGS__nemo__prtctl.flags dom_oce.o in_out_manager.o lib_mpp.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__restart.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

restart.done: restart.o diurnal_bulk.done dom_oce.done eosbn2.done in_out_manager.done iom.done oce.done phycst.done sbc_ice.done trdmxl_oce.done
	touch $(FCM_DONEDIR)/$@

restart.o: $(FCM_PPSRCDIR)/nemo/restart.f90 FFLAGS__nemo__restart.flags diurnal_bulk.o dom_oce.o eosbn2.o in_out_manager.o iom.o oce.o phycst.o sbc_ice.o trdmxl_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbc_ice.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbc_ice.done: sbc_ice.o ice.done in_out_manager.done lib_mpp.done par_oce.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

sbc_ice.o: $(FCM_PPSRCDIR)/nemo/sbc_ice.f90 FFLAGS__nemo__sbc_ice.flags ice.o in_out_manager.o lib_mpp.o par_oce.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbc_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbc_oce.done: sbc_oce.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

sbc_oce.o: $(FCM_PPSRCDIR)/nemo/sbc_oce.f90 FFLAGS__nemo__sbc_oce.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcapr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcapr.done: sbcapr.o dom_oce.done fldread.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done phycst.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

sbcapr.o: $(FCM_PPSRCDIR)/nemo/sbcapr.f90 FFLAGS__nemo__sbcapr.flags dom_oce.o fldread.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o phycst.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcblk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcblk.done: sbcblk.o cyclone.done dom_oce.done fldread.done ice.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limthd_dh.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcblk_algo_coare.done sbcblk_algo_coare3p5.done sbcblk_algo_ecmwf.done sbcblk_algo_ncar.done sbcdcy.done sbcwave.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcblk.o: $(FCM_PPSRCDIR)/nemo/sbcblk.f90 FFLAGS__nemo__sbcblk.flags cyclone.o dom_oce.o fldread.o ice.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limthd_dh.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcblk_algo_coare.o sbcblk_algo_coare3p5.o sbcblk_algo_ecmwf.o sbcblk_algo_ncar.o sbcdcy.o sbcwave.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcblk_algo_coare.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcblk_algo_coare.done: sbcblk_algo_coare.o dom_oce.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcwave.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcblk_algo_coare.o: $(FCM_PPSRCDIR)/nemo/sbcblk_algo_coare.f90 FFLAGS__nemo__sbcblk_algo_coare.flags dom_oce.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcwave.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcblk_algo_coare3p5.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcblk_algo_coare3p5.done: sbcblk_algo_coare3p5.o dom_oce.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcwave.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcblk_algo_coare3p5.o: $(FCM_PPSRCDIR)/nemo/sbcblk_algo_coare3p5.f90 FFLAGS__nemo__sbcblk_algo_coare3p5.flags dom_oce.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcwave.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcblk_algo_ecmwf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcblk_algo_ecmwf.done: sbcblk_algo_ecmwf.o dom_oce.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcwave.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcblk_algo_ecmwf.o: $(FCM_PPSRCDIR)/nemo/sbcblk_algo_ecmwf.f90 FFLAGS__nemo__sbcblk_algo_ecmwf.flags dom_oce.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcwave.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcblk_algo_ncar.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcblk_algo_ncar.done: sbcblk_algo_ncar.o dom_oce.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcwave.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcblk_algo_ncar.o: $(FCM_PPSRCDIR)/nemo/sbcblk_algo_ncar.f90 FFLAGS__nemo__sbcblk_algo_ncar.flags dom_oce.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcwave.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbccpl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbccpl.done: sbccpl.o albedo.done cpl_oasis3.done dom_oce.done eosbn2.done geo2ocean.done ice.done in_out_manager.done iom.done lbclnk.done lib_mpp.done limthd_dh.done oce.done phycst.done sbc_ice.done sbc_oce.done sbcapr.done sbcdcy.done sbcisf.done sbcrnf.done sbcwave.done timing.done trc_oce.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

sbccpl.o: $(FCM_PPSRCDIR)/nemo/sbccpl.f90 FFLAGS__nemo__sbccpl.flags albedo.o cpl_oasis3.o dom_oce.o eosbn2.o geo2ocean.o ice.o in_out_manager.o iom.o lbclnk.o lib_mpp.o limthd_dh.o oce.o phycst.o sbc_ice.o sbc_oce.o sbcapr.o sbcdcy.o sbcisf.o sbcrnf.o sbcwave.o timing.o trc_oce.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcdcy.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcdcy.done: sbcdcy.o dom_oce.done in_out_manager.done lib_mpp.done oce.done phycst.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

sbcdcy.o: $(FCM_PPSRCDIR)/nemo/sbcdcy.f90 FFLAGS__nemo__sbcdcy.flags dom_oce.o in_out_manager.o lib_mpp.o oce.o phycst.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcflx.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcflx.done: sbcflx.o dom_oce.done fldread.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done sbc_oce.done sbcdcy.done
	touch $(FCM_DONEDIR)/$@

sbcflx.o: $(FCM_PPSRCDIR)/nemo/sbcflx.f90 FFLAGS__nemo__sbcflx.flags dom_oce.o fldread.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o sbc_oce.o sbcdcy.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcfwb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcfwb.done: sbcfwb.o dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done sbcisf.done sbcrnf.done sbcssr.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcfwb.o: $(FCM_PPSRCDIR)/nemo/sbcfwb.f90 FFLAGS__nemo__sbcfwb.flags dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o sbcisf.o sbcrnf.o sbcssr.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcice_cice.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcice_cice.done: sbcice_cice.o
	touch $(FCM_DONEDIR)/$@

sbcice_cice.o: $(FCM_PPSRCDIR)/nemo/sbcice_cice.f90 FFLAGS__nemo__sbcice_cice.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcice_if.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcice_if.done: sbcice_if.o dom_oce.done eosbn2.done fldread.done ice.done in_out_manager.done iom.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done
	touch $(FCM_DONEDIR)/$@

sbcice_if.o: $(FCM_PPSRCDIR)/nemo/sbcice_if.f90 FFLAGS__nemo__sbcice_if.flags dom_oce.o eosbn2.o fldread.o ice.o in_out_manager.o iom.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcice_lim.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcice_lim.done: sbcice_lim.o albedo.done bdy_oce.done bdyice_lim.done c1d.done dom_oce.done eosbn2.done ice.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done limctl.done limdiahsb.done limdyn.done limhdf.done limistate.done limitd_me.done limrst.done limsbc.done limthd.done limthd_sal.done limtrp.done limupdate1.done limupdate2.done limvar.done limwri.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcblk.done sbccpl.done thd_ice.done timing.done usrdef_sbc.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcice_lim.o: $(FCM_PPSRCDIR)/nemo/sbcice_lim.f90 FFLAGS__nemo__sbcice_lim.flags albedo.o bdy_oce.o bdyice_lim.o c1d.o dom_oce.o eosbn2.o ice.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o limctl.o limdiahsb.o limdyn.o limhdf.o limistate.o limitd_me.o limrst.o limsbc.o limthd.o limthd_sal.o limtrp.o limupdate1.o limupdate2.o limvar.o limwri.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcblk.o sbccpl.o thd_ice.o timing.o usrdef_sbc.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcice_lim_2.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcice_lim_2.done: sbcice_lim_2.o
	touch $(FCM_DONEDIR)/$@

sbcice_lim_2.o: $(FCM_PPSRCDIR)/nemo/sbcice_lim_2.f90 FFLAGS__nemo__sbcice_lim_2.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcisf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcisf.done: sbcisf.o dom_oce.done eosbn2.done fldread.done in_out_manager.done iom.done lbclnk.done lib_fortran.done oce.done phycst.done sbc_oce.done timing.done wrk_nemo.done zdfbfr.done
	touch $(FCM_DONEDIR)/$@

sbcisf.o: $(FCM_PPSRCDIR)/nemo/sbcisf.f90 FFLAGS__nemo__sbcisf.flags dom_oce.o eosbn2.o fldread.o in_out_manager.o iom.o lbclnk.o lib_fortran.o oce.o phycst.o sbc_oce.o timing.o wrk_nemo.o zdfbfr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcmod.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcmod.done: sbcmod.o bdy_oce.done cpl_oasis3.done diurnal_bulk.done dom_oce.done icbstp.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done prtctl.done sbc_ice.done sbc_oce.done sbcblk.done sbccpl.done sbcdcy.done sbcflx.done sbcfwb.done sbcice_cice.done sbcice_if.done sbcice_lim.done sbcice_lim_2.done sbcisf.done sbcrnf.done sbcssm.done sbcssr.done sbcwave.done timing.done traqsr.done trc_oce.done usrdef_closea.done usrdef_sbc.done
	touch $(FCM_DONEDIR)/$@

sbcmod.o: $(FCM_PPSRCDIR)/nemo/sbcmod.f90 FFLAGS__nemo__sbcmod.flags bdy_oce.o cpl_oasis3.o diurnal_bulk.o dom_oce.o icbstp.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o prtctl.o sbc_ice.o sbc_oce.o sbcblk.o sbccpl.o sbcdcy.o sbcflx.o sbcfwb.o sbcice_cice.o sbcice_if.o sbcice_lim.o sbcice_lim_2.o sbcisf.o sbcrnf.o sbcssm.o sbcssr.o sbcwave.o timing.o traqsr.o trc_oce.o usrdef_closea.o usrdef_sbc.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcrnf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcrnf.done: sbcrnf.o dom_oce.done eosbn2.done fldread.done in_out_manager.done iom.done lib_mpp.done phycst.done sbc_oce.done sbcisf.done usrdef_closea.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sbcrnf.o: $(FCM_PPSRCDIR)/nemo/sbcrnf.f90 FFLAGS__nemo__sbcrnf.flags dom_oce.o eosbn2.o fldread.o in_out_manager.o iom.o lib_mpp.o phycst.o sbc_oce.o sbcisf.o usrdef_closea.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcssm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcssm.done: sbcssm.o c1d.done dom_oce.done eosbn2.done fldread.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done timing.done zdf_oce.done zpshde.done
	touch $(FCM_DONEDIR)/$@

sbcssm.o: $(FCM_PPSRCDIR)/nemo/sbcssm.f90 FFLAGS__nemo__sbcssm.flags c1d.o dom_oce.o eosbn2.o fldread.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o timing.o zdf_oce.o zpshde.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcssr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcssr.done: sbcssr.o dom_oce.done fldread.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done sbcrnf.done timing.done
	touch $(FCM_DONEDIR)/$@

sbcssr.o: $(FCM_PPSRCDIR)/nemo/sbcssr.f90 FFLAGS__nemo__sbcssr.flags dom_oce.o fldread.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o sbcrnf.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbctide.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbctide.done: sbctide.o daymod.done dom_oce.done in_out_manager.done ioipsl.done iom.done lbclnk.done oce.done phycst.done tideini.done
	touch $(FCM_DONEDIR)/$@

sbctide.o: $(FCM_PPSRCDIR)/nemo/sbctide.f90 FFLAGS__nemo__sbctide.flags daymod.o dom_oce.o in_out_manager.o ioipsl.o iom.o lbclnk.o oce.o phycst.o tideini.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sbcwave.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sbcwave.done: sbcwave.o bdy_oce.done domvvl.done fldread.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

sbcwave.o: $(FCM_PPSRCDIR)/nemo/sbcwave.f90 FFLAGS__nemo__sbcwave.flags bdy_oce.o domvvl.o fldread.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__solfrac_mod.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

solfrac_mod.done: solfrac_mod.o par_kind.done
	touch $(FCM_DONEDIR)/$@

solfrac_mod.o: $(FCM_PPSRCDIR)/nemo/solfrac_mod.f90 FFLAGS__nemo__solfrac_mod.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__sshwzv.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

sshwzv.done: sshwzv.o bdy_oce.done bdydyn2d.done divhor.done dom_oce.done domvvl.done in_out_manager.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done restart.done sbc_oce.done timing.done wet_dry.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

sshwzv.o: $(FCM_PPSRCDIR)/nemo/sshwzv.f90 FFLAGS__nemo__sshwzv.flags bdy_oce.o bdydyn2d.o divhor.o dom_oce.o domvvl.o in_out_manager.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o restart.o sbc_oce.o timing.o wet_dry.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__step.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

step.done: step.o bdy_oce.done bdydta.done daymod.done diawri.done dom_oce.done eosbn2.done in_out_manager.done iom.done lbclnk.done oce.done prtctl.done sbc_oce.done sbccpl.done sbcmod.done sbcrnf.done stpctl.done timing.done
	touch $(FCM_DONEDIR)/$@

step.o: $(FCM_PPSRCDIR)/nemo/step.f90 FFLAGS__nemo__step.flags bdy_oce.o bdydta.o daymod.o diawri.o dom_oce.o eosbn2.o in_out_manager.o iom.o lbclnk.o oce.o prtctl.o sbc_oce.o sbccpl.o sbcmod.o sbcrnf.o stpctl.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__step_c1d.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

step_c1d.done: step_c1d.o
	touch $(FCM_DONEDIR)/$@

step_c1d.o: $(FCM_PPSRCDIR)/nemo/step_c1d.f90 FFLAGS__nemo__step_c1d.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__step_diu.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

step_diu.done: step_diu.o cool_skin.done daymod.done diaobs.done diurnal_bulk.done iom.done oce.done restart.done sbc_oce.done sbcmod.done timing.done
	touch $(FCM_DONEDIR)/$@

step_diu.o: $(FCM_PPSRCDIR)/nemo/step_diu.f90 FFLAGS__nemo__step_diu.flags cool_skin.o daymod.o diaobs.o diurnal_bulk.o iom.o oce.o restart.o sbc_oce.o sbcmod.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__step_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

step_oce.done: step_oce.o asmbkg.done asminc.done bdy_oce.done bdydta.done bdydyn3d.done bdytra.done cool_skin.done crsfld.done daymod.done diaar5.done diacfl.done diadct.done diaharm.done diahsb.done diahth.done diaobs.done diaptr.done diawri.done diurnal_bulk.done divhor.done dom_oce.done domvvl.done dynadv.done dynbfr.done dynhpg.done dynldf.done dynnxt.done dynspg.done dynvor.done dynzdf.done eosbn2.done flo_oce.done floats.done in_out_manager.done iom.done lbclnk.done ldfdyn.done ldfslp.done ldftra.done oce.done prtctl.done restart.done sbc_oce.done sbcapr.done sbccpl.done sbcmod.done sbcrnf.done sbctide.done sbcwave.done sshwzv.done step_diu.done stopar.done stopts.done stpctl.done timing.done traadv.done trabbc.done trabbl.done tradmp.done traldf.done tranpc.done tranxt.done traqsr.done trasbc.done trazdf.done zdf_oce.done zdfbfr.done zdfddm.done zdfevd.done zdfgls.done zdfmxl.done zdfqiao.done zdfric.done zdftke.done zdftmx.done zpshde.done
	touch $(FCM_DONEDIR)/$@

step_oce.o: $(FCM_PPSRCDIR)/nemo/step_oce.f90 FFLAGS__nemo__step_oce.flags asmbkg.o asminc.o bdy_oce.o bdydta.o bdydyn3d.o bdytra.o cool_skin.o crsfld.o daymod.o diaar5.o diacfl.o diadct.o diaharm.o diahsb.o diahth.o diaobs.o diaptr.o diawri.o diurnal_bulk.o divhor.o dom_oce.o domvvl.o dynadv.o dynbfr.o dynhpg.o dynldf.o dynnxt.o dynspg.o dynvor.o dynzdf.o eosbn2.o flo_oce.o floats.o in_out_manager.o iom.o lbclnk.o ldfdyn.o ldfslp.o ldftra.o oce.o prtctl.o restart.o sbc_oce.o sbcapr.o sbccpl.o sbcmod.o sbcrnf.o sbctide.o sbcwave.o sshwzv.o step_diu.o stopar.o stopts.o stpctl.o timing.o traadv.o trabbc.o trabbl.o tradmp.o traldf.o tranpc.o tranxt.o traqsr.o trasbc.o trazdf.o zdf_oce.o zdfbfr.o zdfddm.o zdfevd.o zdfgls.o zdfmxl.o zdfqiao.o zdfric.o zdftke.o zdftmx.o zpshde.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__stopar.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

stopar.done: stopar.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done par_oce.done storng.done
	touch $(FCM_DONEDIR)/$@

stopar.o: $(FCM_PPSRCDIR)/nemo/stopar.f90 FFLAGS__nemo__stopar.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o par_oce.o storng.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__stopts.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

stopts.done: stopts.o dom_oce.done lbclnk.done phycst.done stopar.done
	touch $(FCM_DONEDIR)/$@

stopts.o: $(FCM_PPSRCDIR)/nemo/stopts.f90 FFLAGS__nemo__stopts.flags dom_oce.o lbclnk.o phycst.o stopar.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__storng.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

storng.done: storng.o lib_mpp.done par_kind.done
	touch $(FCM_DONEDIR)/$@

storng.o: $(FCM_PPSRCDIR)/nemo/storng.f90 FFLAGS__nemo__storng.flags lib_mpp.o par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__stpctl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

stpctl.done: stpctl.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

stpctl.o: $(FCM_PPSRCDIR)/nemo/stpctl.f90 FFLAGS__nemo__stpctl.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

str_c_to_for.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/str_c_to_for.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

str_c_to_for.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/str_c_to_for.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__thd_ice.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

thd_ice.done: thd_ice.o ice.done in_out_manager.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

thd_ice.o: $(FCM_PPSRCDIR)/nemo/thd_ice.f90 FFLAGS__nemo__thd_ice.flags ice.o in_out_manager.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

tide.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/tide.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

tide.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/tide.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__tide_mod.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

tide_mod.done: tide_mod.o daymod.done dom_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

tide_mod.o: $(FCM_PPSRCDIR)/nemo/tide_mod.f90 FFLAGS__nemo__tide_mod.flags daymod.o dom_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__tideini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

tideini.done: tideini.o daymod.done dom_oce.done in_out_manager.done ioipsl.done iom.done lbclnk.done oce.done phycst.done tide_mod.done
	touch $(FCM_DONEDIR)/$@

tideini.o: $(FCM_PPSRCDIR)/nemo/tideini.f90 FFLAGS__nemo__tideini.flags daymod.o dom_oce.o in_out_manager.o ioipsl.o iom.o lbclnk.o oce.o phycst.o tide_mod.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__timing.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

timing.done: timing.o dom_oce.done in_out_manager.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

timing.o: $(FCM_PPSRCDIR)/nemo/timing.f90 FFLAGS__nemo__timing.flags dom_oce.o in_out_manager.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv.done: traadv.o diaptr.done dom_oce.done domvvl.done in_out_manager.done iom.done ldfslp.done ldftra.done lib_mpp.done oce.done prtctl.done sbc_oce.done sbcwave.done timing.done traadv_cen.done traadv_fct.done traadv_mle.done traadv_mus.done traadv_qck.done traadv_ubs.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv.o: $(FCM_PPSRCDIR)/nemo/traadv.f90 FFLAGS__nemo__traadv.flags diaptr.o dom_oce.o domvvl.o in_out_manager.o iom.o ldfslp.o ldftra.o lib_mpp.o oce.o prtctl.o sbc_oce.o sbcwave.o timing.o traadv_cen.o traadv_fct.o traadv_mle.o traadv_mus.o traadv_qck.o traadv_ubs.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_cen.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_cen.done: traadv_cen.o diaar5.done diaptr.done dom_oce.done eosbn2.done in_out_manager.done iom.done lib_mpp.done oce.done timing.done traadv_fct.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv_cen.o: $(FCM_PPSRCDIR)/nemo/traadv_cen.f90 FFLAGS__nemo__traadv_cen.flags diaar5.o diaptr.o dom_oce.o eosbn2.o in_out_manager.o iom.o lib_mpp.o oce.o timing.o traadv_fct.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_fct.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_fct.done: traadv_fct.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done timing.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv_fct.o: $(FCM_PPSRCDIR)/nemo/traadv_fct.f90 FFLAGS__nemo__traadv_fct.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o timing.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_mle.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_mle.done: traadv_mle.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

traadv_mle.o: $(FCM_PPSRCDIR)/nemo/traadv_mle.f90 FFLAGS__nemo__traadv_mle.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_mus.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_mus.done: traadv_mus.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done oce.done sbcrnf.done timing.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv_mus.o: $(FCM_PPSRCDIR)/nemo/traadv_mus.f90 FFLAGS__nemo__traadv_mus.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o oce.o sbcrnf.o timing.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_qck.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_qck.done: traadv_qck.o diaptr.done dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done oce.done timing.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv_qck.o: $(FCM_PPSRCDIR)/nemo/traadv_qck.f90 FFLAGS__nemo__traadv_qck.flags diaptr.o dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o oce.o timing.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traadv_ubs.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traadv_ubs.done: traadv_ubs.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_fortran.done lib_mpp.done oce.done timing.done traadv_fct.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traadv_ubs.o: $(FCM_PPSRCDIR)/nemo/traadv_ubs.f90 FFLAGS__nemo__traadv_ubs.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_fortran.o lib_mpp.o oce.o timing.o traadv_fct.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trabbc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trabbc.done: trabbc.o dom_oce.done fldread.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done timing.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

trabbc.o: $(FCM_PPSRCDIR)/nemo/trabbc.f90 FFLAGS__nemo__trabbc.flags dom_oce.o fldread.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o timing.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trabbl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trabbl.done: trabbl.o
	touch $(FCM_DONEDIR)/$@

trabbl.o: $(FCM_PPSRCDIR)/nemo/trabbl.f90 FFLAGS__nemo__trabbl.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__tradmp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

tradmp.done: tradmp.o c1d.done dom_oce.done dtatsd.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done prtctl.done timing.done trd_oce.done trdtra.done wrk_nemo.done zdf_oce.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

tradmp.o: $(FCM_PPSRCDIR)/nemo/tradmp.f90 FFLAGS__nemo__tradmp.flags c1d.o dom_oce.o dtatsd.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o prtctl.o timing.o trd_oce.o trdtra.o wrk_nemo.o zdf_oce.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traldf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traldf.done: traldf.o dom_oce.done in_out_manager.done lbclnk.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done prtctl.done timing.done traldf_iso.done traldf_lap_blp.done traldf_triad.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traldf.o: $(FCM_PPSRCDIR)/nemo/traldf.f90 FFLAGS__nemo__traldf.flags dom_oce.o in_out_manager.o lbclnk.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o prtctl.o timing.o traldf_iso.o traldf_lap_blp.o traldf_triad.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traldf_iso.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traldf_iso.done: traldf_iso.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done ldfslp.done ldftra.done oce.done phycst.done timing.done trc_oce.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

traldf_iso.o: $(FCM_PPSRCDIR)/nemo/traldf_iso.f90 FFLAGS__nemo__traldf_iso.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o ldfslp.o ldftra.o oce.o phycst.o timing.o trc_oce.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traldf_lap_blp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traldf_lap_blp.done: traldf_lap_blp.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done ldftra.done lib_mpp.done oce.done timing.done traldf_iso.done traldf_triad.done trc_oce.done wrk_nemo.done zpshde.done
	touch $(FCM_DONEDIR)/$@

traldf_lap_blp.o: $(FCM_PPSRCDIR)/nemo/traldf_lap_blp.f90 FFLAGS__nemo__traldf_lap_blp.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o ldftra.o lib_mpp.o oce.o timing.o traldf_iso.o traldf_triad.o trc_oce.o wrk_nemo.o zpshde.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traldf_triad.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traldf_triad.done: traldf_triad.o diaar5.done diaptr.done dom_oce.done in_out_manager.done iom.done lbclnk.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done timing.done traldf_iso.done trc_oce.done wrk_nemo.done zdf_oce.done zpshde.done
	touch $(FCM_DONEDIR)/$@

traldf_triad.o: $(FCM_PPSRCDIR)/nemo/traldf_triad.f90 FFLAGS__nemo__traldf_triad.flags diaar5.o diaptr.o dom_oce.o in_out_manager.o iom.o lbclnk.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o timing.o traldf_iso.o trc_oce.o wrk_nemo.o zdf_oce.o zpshde.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__tranpc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

tranpc.done: tranpc.o dom_oce.done eosbn2.done in_out_manager.done lbclnk.done lib_mpp.done oce.done phycst.done timing.done trd_oce.done trdtra.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

tranpc.o: $(FCM_PPSRCDIR)/nemo/tranpc.f90 FFLAGS__nemo__tranpc.flags dom_oce.o eosbn2.o in_out_manager.o lbclnk.o lib_mpp.o oce.o phycst.o timing.o trd_oce.o trdtra.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__tranxt.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

tranxt.done: tranxt.o bdy_oce.done bdytra.done dom_oce.done domvvl.done in_out_manager.done lbclnk.done ldfslp.done ldftra.done oce.done phycst.done prtctl.done sbc_oce.done sbcisf.done sbcrnf.done timing.done traqsr.done trd_oce.done trdtra.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

tranxt.o: $(FCM_PPSRCDIR)/nemo/tranxt.f90 FFLAGS__nemo__tranxt.flags bdy_oce.o bdytra.o dom_oce.o domvvl.o in_out_manager.o lbclnk.o ldfslp.o ldftra.o oce.o phycst.o prtctl.o sbc_oce.o sbcisf.o sbcrnf.o timing.o traqsr.o trd_oce.o trdtra.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__traqsr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

traqsr.done: traqsr.o dom_oce.done fldread.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done restart.done sbc_oce.done timing.done trc_oce.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

traqsr.o: $(FCM_PPSRCDIR)/nemo/traqsr.f90 FFLAGS__nemo__traqsr.flags dom_oce.o fldread.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o restart.o sbc_oce.o timing.o trc_oce.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trasbc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trasbc.done: trasbc.o dom_oce.done eosbn2.done in_out_manager.done iom.done iscplini.done lbclnk.done oce.done phycst.done prtctl.done sbc_oce.done sbcisf.done sbcmod.done sbcrnf.done timing.done traqsr.done trd_oce.done trdtra.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

trasbc.o: $(FCM_PPSRCDIR)/nemo/trasbc.f90 FFLAGS__nemo__trasbc.flags dom_oce.o eosbn2.o in_out_manager.o iom.o iscplini.o lbclnk.o oce.o phycst.o prtctl.o sbc_oce.o sbcisf.o sbcmod.o sbcrnf.o timing.o traqsr.o trd_oce.o trdtra.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trazdf.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trazdf.done: trazdf.o dom_oce.done domvvl.done in_out_manager.done lbclnk.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done prtctl.done sbc_oce.done timing.done trazdf_exp.done trazdf_imp.done trd_oce.done trdtra.done wrk_nemo.done zdf_oce.done zdfgls.done zdftke.done
	touch $(FCM_DONEDIR)/$@

trazdf.o: $(FCM_PPSRCDIR)/nemo/trazdf.f90 FFLAGS__nemo__trazdf.flags dom_oce.o domvvl.o in_out_manager.o lbclnk.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o prtctl.o sbc_oce.o timing.o trazdf_exp.o trazdf_imp.o trd_oce.o trdtra.o wrk_nemo.o zdf_oce.o zdfgls.o zdftke.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trazdf_exp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trazdf_exp.done: trazdf_exp.o dom_oce.done domvvl.done in_out_manager.done lib_mpp.done oce.done timing.done trc_oce.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

trazdf_exp.o: $(FCM_PPSRCDIR)/nemo/trazdf_exp.f90 FFLAGS__nemo__trazdf_exp.flags dom_oce.o domvvl.o in_out_manager.o lib_mpp.o oce.o timing.o trc_oce.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trazdf_imp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trazdf_imp.done: trazdf_imp.o dom_oce.done domvvl.done in_out_manager.done lbclnk.done ldfslp.done ldftra.done lib_mpp.done oce.done timing.done traldf_triad.done trc_oce.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

trazdf_imp.o: $(FCM_PPSRCDIR)/nemo/trazdf_imp.f90 FFLAGS__nemo__trazdf_imp.flags dom_oce.o domvvl.o in_out_manager.o lbclnk.o ldfslp.o ldftra.o lib_mpp.o oce.o timing.o traldf_triad.o trc_oce.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trc_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trc_oce.done: trc_oce.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

trc_oce.o: $(FCM_PPSRCDIR)/nemo/trc_oce.f90 FFLAGS__nemo__trc_oce.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trd_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trd_oce.done: trd_oce.o par_oce.done trdmxl_oce.done trdvor_oce.done
	touch $(FCM_DONEDIR)/$@

trd_oce.o: $(FCM_PPSRCDIR)/nemo/trd_oce.f90 FFLAGS__nemo__trd_oce.flags par_oce.o trdmxl_oce.o trdvor_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trddyn.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trddyn.done: trddyn.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done trdglo.done trdken.done trdmxl.done trdvor.done wrk_nemo.done zdf_oce.done zdfbfr.done
	touch $(FCM_DONEDIR)/$@

trddyn.o: $(FCM_PPSRCDIR)/nemo/trddyn.f90 FFLAGS__nemo__trddyn.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o trdglo.o trdken.o trdmxl.o trdvor.o wrk_nemo.o zdf_oce.o zdfbfr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdglo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdglo.done: trdglo.o dom_oce.done eosbn2.done in_out_manager.done iom.done ldfdyn.done ldftra.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done wrk_nemo.done zdf_oce.done zdfbfr.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

trdglo.o: $(FCM_PPSRCDIR)/nemo/trdglo.f90 FFLAGS__nemo__trdglo.flags dom_oce.o eosbn2.o in_out_manager.o iom.o ldfdyn.o ldftra.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o wrk_nemo.o zdf_oce.o zdfbfr.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdini.done: trdini.o dom_oce.done in_out_manager.done lib_mpp.done trd_oce.done trdglo.done trdken.done trdmxl.done trdpen.done trdvor.done
	touch $(FCM_DONEDIR)/$@

trdini.o: $(FCM_PPSRCDIR)/nemo/trdini.f90 FFLAGS__nemo__trdini.flags dom_oce.o in_out_manager.o lib_mpp.o trd_oce.o trdglo.o trdken.o trdmxl.o trdpen.o trdvor.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdken.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdken.done: trdken.o dom_oce.done in_out_manager.done iom.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done trdglo.done trdmxl.done trdvor.done wrk_nemo.done zdf_oce.done zdfbfr.done
	touch $(FCM_DONEDIR)/$@

trdken.o: $(FCM_PPSRCDIR)/nemo/trdken.f90 FFLAGS__nemo__trdken.flags dom_oce.o in_out_manager.o iom.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o trdglo.o trdmxl.o trdvor.o wrk_nemo.o zdf_oce.o zdfbfr.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdmxl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdmxl.done: trdmxl.o dianam.done dom_oce.done in_out_manager.done ioipsl.done iom.done lbclnk.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done prtctl.done restart.done trd_oce.done trdmxl_oce.done trdmxl_rst.done wrk_nemo.done zdf_oce.done zdfddm.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

trdmxl.o: $(FCM_PPSRCDIR)/nemo/trdmxl.f90 FFLAGS__nemo__trdmxl.flags dianam.o dom_oce.o in_out_manager.o ioipsl.o iom.o lbclnk.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o prtctl.o restart.o trd_oce.o trdmxl_oce.o trdmxl_rst.o wrk_nemo.o zdf_oce.o zdfddm.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdmxl_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdmxl_oce.done: trdmxl_oce.o lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

trdmxl_oce.o: $(FCM_PPSRCDIR)/nemo/trdmxl_oce.f90 FFLAGS__nemo__trdmxl_oce.flags lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdmxl_rst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdmxl_rst.done: trdmxl_rst.o dom_oce.done in_out_manager.done iom.done restart.done trd_oce.done
	touch $(FCM_DONEDIR)/$@

trdmxl_rst.o: $(FCM_PPSRCDIR)/nemo/trdmxl_rst.f90 FFLAGS__nemo__trdmxl_rst.flags dom_oce.o in_out_manager.o iom.o restart.o trd_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdpen.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdpen.done: trdpen.o dom_oce.done eosbn2.done in_out_manager.done iom.done ldftra.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

trdpen.o: $(FCM_PPSRCDIR)/nemo/trdpen.f90 FFLAGS__nemo__trdpen.flags dom_oce.o eosbn2.o in_out_manager.o iom.o ldftra.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdtra.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdtra.done: trdtra.o dom_oce.done in_out_manager.done iom.done ldfslp.done ldftra.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done trdglo.done trdmxl.done trdpen.done trdtrc.done wrk_nemo.done zdf_oce.done zdfddm.done
	touch $(FCM_DONEDIR)/$@

trdtra.o: $(FCM_PPSRCDIR)/nemo/trdtra.f90 FFLAGS__nemo__trdtra.flags dom_oce.o in_out_manager.o iom.o ldfslp.o ldftra.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o trdglo.o trdmxl.o trdpen.o trdtrc.o wrk_nemo.o zdf_oce.o zdfddm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdtrc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdtrc.done: trdtrc.o
	touch $(FCM_DONEDIR)/$@

trdtrc.o: $(FCM_PPSRCDIR)/nemo/trdtrc.f90 FFLAGS__nemo__trdtrc.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdvor.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdvor.done: trdvor.o dianam.done dom_oce.done in_out_manager.done ioipsl.done lbclnk.done ldfdyn.done lib_mpp.done oce.done phycst.done sbc_oce.done trd_oce.done wrk_nemo.done zdf_oce.done zdfmxl.done
	touch $(FCM_DONEDIR)/$@

trdvor.o: $(FCM_PPSRCDIR)/nemo/trdvor.f90 FFLAGS__nemo__trdvor.flags dianam.o dom_oce.o in_out_manager.o ioipsl.o lbclnk.o ldfdyn.o lib_mpp.o oce.o phycst.o sbc_oce.o trd_oce.o wrk_nemo.o zdf_oce.o zdfmxl.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trdvor_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trdvor_oce.done: trdvor_oce.o par_oce.done
	touch $(FCM_DONEDIR)/$@

trdvor_oce.o: $(FCM_PPSRCDIR)/nemo/trdvor_oce.f90 FFLAGS__nemo__trdvor_oce.flags par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__updtide.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

updtide.done: updtide.o dom_oce.done in_out_manager.done oce.done phycst.done sbctide.done tideini.done
	touch $(FCM_DONEDIR)/$@

updtide.o: $(FCM_PPSRCDIR)/nemo/updtide.f90 FFLAGS__nemo__updtide.flags dom_oce.o in_out_manager.o oce.o phycst.o sbctide.o tideini.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_closea.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_closea.done: usrdef_closea.o dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done timing.done
	touch $(FCM_DONEDIR)/$@

usrdef_closea.o: $(FCM_PPSRCDIR)/nemo/usrdef_closea.f90 FFLAGS__nemo__usrdef_closea.flags dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_fmask.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_fmask.done: usrdef_fmask.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

usrdef_fmask.o: $(FCM_PPSRCDIR)/nemo/usrdef_fmask.f90 FFLAGS__nemo__usrdef_fmask.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_hgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_hgr.done: usrdef_hgr.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done phycst.done timing.done usrdef_nam.done
	touch $(FCM_DONEDIR)/$@

usrdef_hgr.o: $(FCM_PPSRCDIR)/nemo/usrdef_hgr.f90 FFLAGS__nemo__usrdef_hgr.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o phycst.o timing.o usrdef_nam.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_istate.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_istate.done: usrdef_istate.o in_out_manager.done lib_mpp.done par_oce.done phycst.done
	touch $(FCM_DONEDIR)/$@

usrdef_istate.o: $(FCM_PPSRCDIR)/nemo/usrdef_istate.f90 FFLAGS__nemo__usrdef_istate.flags in_out_manager.o lib_mpp.o par_oce.o phycst.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_nam.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_nam.done: usrdef_nam.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

usrdef_nam.o: $(FCM_PPSRCDIR)/nemo/usrdef_nam.f90 FFLAGS__nemo__usrdef_nam.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_sbc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_sbc.done: usrdef_sbc.o dom_oce.done ice.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done limthd_dh.done oce.done phycst.done sbc_ice.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

usrdef_sbc.o: $(FCM_PPSRCDIR)/nemo/usrdef_sbc.f90 FFLAGS__nemo__usrdef_sbc.flags dom_oce.o ice.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o limthd_dh.o oce.o phycst.o sbc_ice.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__usrdef_zgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

usrdef_zgr.done: usrdef_zgr.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done timing.done usrdef_nam.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

usrdef_zgr.o: $(FCM_PPSRCDIR)/nemo/usrdef_zgr.f90 FFLAGS__nemo__usrdef_zgr.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o timing.o usrdef_nam.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

vectopt_loop_substitute.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/vectopt_loop_substitute.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

vectopt_loop_substitute.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/vectopt_loop_substitute.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__wet_dry.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

wet_dry.done: wet_dry.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done sbc_oce.done sbcrnf.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

wet_dry.o: $(FCM_PPSRCDIR)/nemo/wet_dry.f90 FFLAGS__nemo__wet_dry.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o sbc_oce.o sbcrnf.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__wrk_nemo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

wrk_nemo.done: wrk_nemo.o par_oce.done
	touch $(FCM_DONEDIR)/$@

wrk_nemo.o: $(FCM_PPSRCDIR)/nemo/wrk_nemo.f90 FFLAGS__nemo__wrk_nemo.flags par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdf_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdf_oce.done: zdf_oce.o in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

zdf_oce.o: $(FCM_PPSRCDIR)/nemo/zdf_oce.f90 FFLAGS__nemo__zdf_oce.flags in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfbfr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfbfr.done: zdfbfr.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done prtctl.done timing.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

zdfbfr.o: $(FCM_PPSRCDIR)/nemo/zdfbfr.f90 FFLAGS__nemo__zdfbfr.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o prtctl.o timing.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfddm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfddm.done: zdfddm.o
	touch $(FCM_DONEDIR)/$@

zdfddm.o: $(FCM_PPSRCDIR)/nemo/zdfddm.f90 FFLAGS__nemo__zdfddm.flags
	fcm_internal compile:F nemo $< $@

zdfddm_substitute.h90: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/zdfddm_substitute.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

zdfddm_substitute.h90.idone: /home1/datawork/gboutin/trunk/NEMOGCM/CONFIG/TEST_CASES/clean_with_clement/WORK/zdfddm_substitute.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__zdfevd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfevd.done: zdfevd.o dom_oce.done in_out_manager.done iom.done lbclnk.done oce.done timing.done trd_oce.done trdtra.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

zdfevd.o: $(FCM_PPSRCDIR)/nemo/zdfevd.f90 FFLAGS__nemo__zdfevd.flags dom_oce.o in_out_manager.o iom.o lbclnk.o oce.o timing.o trd_oce.o trdtra.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfgls.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfgls.done: zdfgls.o
	touch $(FCM_DONEDIR)/$@

zdfgls.o: $(FCM_PPSRCDIR)/nemo/zdfgls.f90 FFLAGS__nemo__zdfgls.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfini.done: zdfini.o in_out_manager.done iom.done ldfslp.done lib_mpp.done par_oce.done sbc_oce.done tranpc.done zdf_oce.done zdfddm.done zdfevd.done zdfgls.done zdfric.done zdftke.done
	touch $(FCM_DONEDIR)/$@

zdfini.o: $(FCM_PPSRCDIR)/nemo/zdfini.f90 FFLAGS__nemo__zdfini.flags in_out_manager.o iom.o ldfslp.o lib_mpp.o par_oce.o sbc_oce.o tranpc.o zdf_oce.o zdfddm.o zdfevd.o zdfgls.o zdfric.o zdftke.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfmxl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfmxl.done: zdfmxl.o dom_oce.done in_out_manager.done iom.done lib_mpp.done oce.done phycst.done prtctl.done timing.done trc_oce.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

zdfmxl.o: $(FCM_PPSRCDIR)/nemo/zdfmxl.f90 FFLAGS__nemo__zdfmxl.flags dom_oce.o in_out_manager.o iom.o lib_mpp.o oce.o phycst.o prtctl.o timing.o trc_oce.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfqiao.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfqiao.done: zdfqiao.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done sbc_oce.done sbcwave.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

zdfqiao.o: $(FCM_PPSRCDIR)/nemo/zdfqiao.f90 FFLAGS__nemo__zdfqiao.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o sbc_oce.o sbcwave.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdfric.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdfric.done: zdfric.o
	touch $(FCM_DONEDIR)/$@

zdfric.o: $(FCM_PPSRCDIR)/nemo/zdfric.f90 FFLAGS__nemo__zdfric.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdftke.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdftke.done: zdftke.o
	touch $(FCM_DONEDIR)/$@

zdftke.o: $(FCM_PPSRCDIR)/nemo/zdftke.f90 FFLAGS__nemo__zdftke.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zdftmx.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zdftmx.done: zdftmx.o
	touch $(FCM_DONEDIR)/$@

zdftmx.o: $(FCM_PPSRCDIR)/nemo/zdftmx.f90 FFLAGS__nemo__zdftmx.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zpshde.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zpshde.done: zpshde.o dom_oce.done eosbn2.done in_out_manager.done lbclnk.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

zpshde.o: $(FCM_PPSRCDIR)/nemo/zpshde.f90 FFLAGS__nemo__zpshde.flags dom_oce.o eosbn2.o in_out_manager.o lbclnk.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

# EOF
