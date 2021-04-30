// **************************************************************************
//                               charmm.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for acceleration of the charmm/coul pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture(pos_tex, float4);
_texture(q_tex, float);
#else
_texture_2d(pos_tex, int4);
_texture(q_tex, int2);
#endif

#else
#define pos_tex x_
#define q_tex q_
#endif

__kernel void k_charmm(const __global numtyp4 *restrict x_,
                       const __global numtyp4 *restrict lj1,
                       const int lj_types,
                       const __global numtyp *restrict sp_lj,
                       const __global int *dev_nbor,
                       const __global int *dev_packed,
                       __global acctyp4 *restrict ans,
                       __global acctyp *restrict engv,
                       const int eflag, const int vflag,
                       const int inum, const int nbor_pitch,
                       const __global numtyp *restrict q_,
                       const numtyp cut_coulsq, const numtyp qqrd2e,
                       const numtyp denom_lj,
                       const numtyp denom_coul,
                       const numtyp cut_bothsq,
                       const numtyp cut_ljsq,
                       const numtyp cut_lj_innersq,
                       const numtyp cut_coul_innersq,
                       const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_bio();

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cut_bothsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force_lj, force, r6inv, switch1;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            numtyp switch2 = (numtyp)12.0*rsq*switch1*(rsq-cut_lj_innersq)*
                             denom_lj;
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)*
                       denom_lj;
            switch2 *= r6inv*(lj1[mtype].z*r6inv-lj1[mtype].w);
            force_lj = force_lj*switch1+switch2;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp rinv = ucl_rsqrt(rsq);
          fetch(forcecoul,j,q_tex);
          forcecoul *= factor_coul * qqrd2e * qtmp * rinv;
          if (rsq > cut_coul_innersq) {
            numtyp switch3 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
              (cut_coulsq + (numtyp)2.0*rsq - (numtyp)3.0*cut_coul_innersq) *
              denom_coul;
            forcecoul *= switch3;
          }
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += forcecoul;
          if (rsq < cut_ljsq) {
            numtyp e=r6inv*(lj1[mtype].z*r6inv-lj1[mtype].w);
            if (rsq > cut_lj_innersq)
              e *= switch1;
            energy+=factor_lj*e;
          }
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

__kernel void k_charmm_fast(const __global numtyp4 *restrict x_,
                            const __global numtyp2 *restrict ljd_in,
                            const __global numtyp *restrict sp_lj_in,
                            const __global int *dev_nbor,
                            const __global int *dev_packed,
                            __global acctyp4 *restrict ans,
                            __global acctyp *restrict engv,
                            const int eflag, const int vflag,
                            const int inum, const int nbor_pitch,
                            const __global numtyp *restrict q_,
                            const numtyp cut_coulsq, const numtyp qqrd2e,
                            const numtyp denom_lj,
                            const numtyp denom_coul,
                            const numtyp cut_bothsq,
                            const numtyp cut_ljsq,
                            const numtyp cut_lj_innersq,
                            const numtyp cut_coul_innersq,
                            const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp2 ljd[MAX_BIO_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_bio();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_BIO_SHARED_TYPES)
    ljd[tid]=ljd_in[tid];
  if (tid+BLOCK_BIO_PAIR<MAX_BIO_SHARED_TYPES)
    ljd[tid+BLOCK_BIO_PAIR]=ljd_in[tid+BLOCK_BIO_PAIR];

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cut_bothsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force_lj, force, switch1;
        numtyp lj3, lj4;

        if (rsq < cut_ljsq) {
          numtyp eps = ucl_sqrt(ljd[itype].x*ljd[jtype].x);
          numtyp sig6 = (numtyp)0.5 * (ljd[itype].y+ljd[jtype].y);

          numtyp sig_r_6 = sig6*sig6*r2inv;
          sig_r_6 = sig_r_6*sig_r_6*sig_r_6;
          lj4 = (numtyp)4.0*eps*sig_r_6;
          lj3 = lj4*sig_r_6;
          force_lj = factor_lj*((numtyp)12.0 * lj3 - (numtyp)6.0 * lj4);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            numtyp switch2 = (numtyp)12.0*rsq*switch1*(rsq-cut_lj_innersq)*
                             denom_lj;
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)*
                       denom_lj;
            switch2 *= lj3-lj4;
            force_lj = force_lj*switch1+switch2;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp rinv = ucl_rsqrt(rsq);
          fetch(forcecoul,j,q_tex);
          forcecoul *= factor_coul * qqrd2e * qtmp * rinv;
          if (rsq > cut_coul_innersq) {
            numtyp switch3 = (cut_coulsq-rsq) * (cut_coulsq-rsq) *
              (cut_coulsq + (numtyp)2.0*rsq - (numtyp)3.0*cut_coul_innersq) *
              denom_coul;
            forcecoul *= switch3;
          }
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += forcecoul;
          if (rsq < cut_ljsq) {
            numtyp e=lj3-lj4;
            if (rsq > cut_lj_innersq)
              e *= switch1;
            energy+=factor_lj*e;
          }
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}
