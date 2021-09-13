 
     SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
     !------------------------------------------------------------
     !This subroutine is called at the very start of the analysis-
     !and is used to read in the matrices containing all grain   -
     !boundary information.                                      -
     !Replace the location to the files as required.             -     
     !------------------------------------------------------------

      use ISO_FORTRAN_ENV
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'param_array.inc'

      DIMENSION TIME(2)

      
      integer*8 :: nr,nc,intotalfeat
      real*8 :: numrowval,numcolval,
     1 boundgrain(totalfeat,nodeout),elcent(totalels,3),
     2 nodex(totalfeat,nodeout),nodey(totalfeat,nodeout),
     3 nodez(totalfeat,nodeout)

      real*8, allocatable :: rowdata(:)
	
      common nodex,nodey,nodez,boundgrain,elcent,intotalfeat

       if(lop .eq. 0)then

      	! Read in bound feature array
      	open(500, FILE="/.../boundfeat.bin", 		    
     1    form='unformatted',status='old',access='stream')
                
      	read(500) numrowval
      	nr=int(numrowval)
     
      	read(500) numcolval
      	nc=int(numcolval)
     
      	allocate (rowdata(nc))
 
      	do i=1,totalfeat
        	read(500) rowdata
        	boundgrain(i,:)=rowdata(:)
      	end do
	
      	deallocate(rowdata)

      	!Read in element centroid coordinates
      	open(600, FILE="/.../el_centroid.bin", 		    
     1  form='unformatted',status='old',access='stream')
                
                
      	read(600) numrowval
      	nr=int(numrowval)
     
      	read(600) numcolval
      	nc=int(numcolval)
     
      	allocate (rowdata(nc))
     
      	do i=1,totalels
        	read(600) rowdata
        	elcent(i,:)=rowdata(:)
      	end do
     
      
      	deallocate(rowdata)

      	! reading in the x, y, and z coordinates for nodes on the outside
      	! of the grains
       	open(233, FILE="/.../xvalues.bin", 	   
     1  form='unformatted',status='old',access='stream')
                
      	read(233) numrowval
      	nr=int(numrowval)

      	read(233) numcolval
      	nc=int(numcolval)
	
      	open(244, FILE="/.../yvalues.bin", 	   
     1  form='unformatted',status='old',access='stream')

      	read(244) numrowval
      	nr=int(numrowval)
     
      	read(244) numcolval
      	nc=int(numcolval)
                     
      	open(225, FILE="/.../zvalues.bin",	   
     1    form='unformatted', status='old',access='stream')

      	read(225) numrowval
      	nr=int(numrowval)
     
      	read(225) numcolval
      	nc=int(numcolval)
     
      	allocate (rowdata(nc))
     
      
      	do i=1,nr
        	read(233) rowdata
        	nodex(i,:)=rowdata(:)
	
        
        	read(244) rowdata
        	nodey(i,:)=rowdata(:)
        
        	read(225) rowdata
        	nodez(i,:)=rowdata(:)
       
      	end do
        deallocate(rowdata)
      	intotalfeat=nr


		close(233)
		close(244)
		close(225)
		close(600)
		close(500)

       	end if

       RETURN
       end subroutine  UEXTERNALDB
