      program ifdef
c
      character*80 c(2)
      character*3 v
c
      j1=1
      j2=2
c
    1 format(a80)
    2 format(a7,a3)
    3 format(a5,a75)
    4 format(a6)
c
      read(1,1) c(j1)
      if(c(j1)(1:2).eq.'##') then
        v=c(j1)(3:5)
        do i=1,3
          call l2u(v(i:i))
        enddo
        write(51,2) '#ifdef ',v
        write(51,3) '     ',c(j1)(6:80)
      else
        write(51,1) c(j1)
      endif
c
  100 continue
      read(1,1,end=200) c(j2)
      if(c(j1)(1:2).ne.'##') then
        if(c(j2)(1:2).eq.'##') then
          v=c(j2)(3:5)
          do i=1,3
            call l2u(v(i:i))
          enddo
          write(51,2) '#ifdef ',v
          write(51,3) '     ',c(j2)(6:80)
        else
          write(51,1) c(j2)
        endif
      else
        if(c(j2)(1:2).eq.'##') then
          if(c(j2)(1:5).eq.c(j1)(1:5)) then
            write(51,3) '     ',c(j2)(6:80)
          else
            write(51,4) '#endif'
            v=c(j2)(3:5)
            do i=1,3
              call l2u(v(i:i))
            enddo
            write(51,2) '#ifdef ',v
            write(51,3) '     ',c(j2)(6:80)
          endif
        else
          write(51,4) '#endif'
          write(51,1) c(j2)
        endif
      endif
  101 continue
      jx=j1
      j1=j2
      j2=jx
      go to 100
  200 continue
      if(c(j1)(1:2).eq.'##') then
        write(51,4) '#endif'
      endif
      stop
      end
      subroutine l2u(c)
c
      character*1 c
c
      character*1 lower(26),upper(26)
c
      data lower/'a','b','c','d','e','f','g','h','i','j','k',
     1           'l','m','n','o','p','q','r','s','t','u','v',
     2           'w','x','y','z'/
      data upper/'A','B','C','D','E','F','G','H','I','J','K',
     1           'L','M','N','O','P','Q','R','S','T','U','V',
     2           'W','X','Y','Z'/
c
      do i=1,26
        if(c.eq.lower(i)) then
          c=upper(i)
          return
        endif
      enddo
      return
      end
