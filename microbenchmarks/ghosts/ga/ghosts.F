#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#define MPI
      program main

      integer NDIM, NSIZE, ITERS
      parameter (NDIM=2, NSIZE=513, ITERS=5000)
      common /indices/ size(NDIM), width(NDIM), dims(NDIM)
      integer size, width, dims
      common /ghstdat/ dims_fg(NDIM+1)
      integer dims_fg
      integer i, j, gdim, type
      integer gdims(NDIM+1), gwidth(NDIM+1)
      integer pdims(NDIM+1), mcnt, mapc(5000)
      integer g_fg, ptr_fg
      integer ld_fg(NDIM+1), ld_fld(NDIM+1), ld_bc(NDIM)
      integer heap, stack, me, nproc
      double precision start, time
c
#include "mafdecls.fh"
#include "global.fh"
c
c  Initialize a message passing library
c
#include "mp3.fh"
      size(1) = NSIZE
      size(2) = NSIZE
c
c   Initialize global arrays
c
      call ga_initialize()
c
      nproc = ga_nnodes()
      me = ga_nodeid()
c
      if (ga_uses_ma()) then
        heap = (size(1)+2)*(size(2)+2)*34/nproc
      else
        heap = 100000
      endif
      stack = 50000
c
      if (.not.ma_init(MT_DBL, stack, heap))
     +  call ga_error("ma init failed", -1)

c
c initialize global arrays
c
      do i = 1, NDIM
        dims(i) = size(i)
        gdims(i) = size(i)
        width(i) = 1
        gwidth(i) = 1
      end do
c
c  evaluate distribution of processors
c
      gdim = NDIM
      call factor(nproc,gdim,pdims)
      mcnt = 1
      do i = 1, NDIM
        do j = 0, pdims(i)-1
          mapc(mcnt) = ((dble(j)/dble(pdims(i)))*dble(gdims(i)))+1
          mcnt = mcnt + 1
        end do
      end do
      do i = 0, pdims(1)-1
        mapc(mcnt) = ((dble(i)/dble(pdims(1)))*dble(NDIM))+1
        mcnt = mcnt + 1
      end do
      do i = 0, pdims(2)-1
        mapc(mcnt) = ((dble(i)/dble(pdims(2)))*dble(NDIM))+1
        mcnt = mcnt + 1
      end do
c
c   Create global array
c
      type = MT_DBL
      gdim = 3
      gdims(3) = 9
      gwidth(3) = 0
      pdims(3) = 1
      mapc(mcnt) = 1
      if (.not.nga_create_ghosts_irreg(type, gdim, gdims, gwidth,
     +    "lb_dist", mapc, pdims, g_fg))
     +    call ga_error("g_fg init failed",me)
c
c   Find pointers to global array data
c
      call nga_access_ghosts(g_fg,dims_fg,ptr_fg,ld_fg)
      do i = 1, NDIM
        dims(i) = dims_fg(i)
      end do
c
c   Call routine to run main simulation
c
      start = MP_TIMER()
      do i = 1, ITERS
      call ga_mask_sync(.false.,.false.)
      call ga_update_ghosts(g_fg)
      end do
      time = MP_TIMER() - start

      if (me.eq.0) write(6,100) (time*1000d0)/ITERS

  100 format('updateGhosts avg :',f12.4,'ms')
c
c   Close out calculation
c
      call ga_terminate()
      call MP_FINALIZE()
      stop
      end
c
      subroutine factor(p,ndim,dims)
      implicit none
      integer i,j,p,ndim,dims(*),imin,mdim
      integer ip,ifac,pmax,prime(1000)
      integer fac(1000)
c
      i = 1
      ip = p
      do i = 1, ndim
        dims(i) = 1
      end do
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine dimensions of processor grid
c
      do i = ifac, 1, -1
c
c    find dimension with minimum value
c
        imin = dims(1)
        mdim = 1
        do j = 2, ndim
          if (dims(j).lt.imin) then
            imin = dims(j)
            mdim = j
          endif
        end do
        dims(mdim) = dims(mdim)*fac(i)
      end do
c
      return
      end
