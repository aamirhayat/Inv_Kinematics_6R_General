
       subroutine read_a(a1,a2,a3,a4,a5,a6)

       double precision a1, a2, a3, a4, a5, a6

       read(5,*) a1, a2, a3, a4, a5, a6
       write(6,*) a1, a2, a3, a4, a5, a6
      
       return
       end

       subroutine read_d(d1,d2,d3,d4,d5,d6)

       double precision d1, d2, d3, d4, d5, d6

       read(5,*) d1, d2, d3, d4, d5, d6
       write(6,*) d1, d2, d3, d4, d5, d6
       return
       end


       subroutine read_al(al1,al2,al3,al4,al5,al6)

       double precision al1,al2,al3,al4,al5,al6

       read(5,*) al1,al2,al3,al4,al5,al6
       write(6,*) al1,al2,al3,al4,al5,al6
       return
       end


       subroutine read_rhs(lx,mx,nx,ly,my,ny,lz,mz,nz,qx,qy,qz)

       double precision lx,mx,nx,ly,my,ny,lz,mz,nz,qx,qy,qz

       read(5,*)lx,mx,nx,ly,my,ny,lz,mz,nz,qx,qy,qz
       write(6,*)lx,mx,nx,ly,my,ny,lz,mz,nz,qx,qy,qz
       return
       end

