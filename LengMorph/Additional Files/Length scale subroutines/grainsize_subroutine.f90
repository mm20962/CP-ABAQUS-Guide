!----------------------------------------------------------------------------------------
!subroutine containing the length scale calculation which can be added to crystal plasticity
!models.
subroutine lengthscale(coords,nodeout,totalfeat,nodex,nodey,nodez,      &
     slpdir1,slpnor1,nslptl,elcent,feature,boundgrain,oriensh,noel,   &
     rdistance,ldistance,lm)

!INPUTS :
!coords             - voxel coordinates
!nodeout            - columns in nodex, nodey, nodez,boundfeat
!totalfeat          - total number of grains in the RVE
!nodex,nodey,nodez  - x, y, and z coordinates for each boundary node
!slpdir1             - slip directions in local (crystal) coordinate system
!slpnor1             - slip plane normals in local (crystal) coordinate system
!elcent             - array containing the centroids of all voxels
!feature            - grain number
!boundgrain         - boundary grain ids
!oriensh            - orientations in euler angles in each grain
!noel               - voxel number

!OUTPUTS:
!rdistance          - distance from the voxel to grain boundary
!ldistance          - distance along slip length
!lm                 - Luster-Morris parameter


IMPLICIT REAL*8 (A-H,O-Z)
parameter(nd=150)
integer :: i, arraysize, nodeout,totalfeat,feature,tester,grainb,&
    arraysizeb,val,noel,valb
integer*8 :: minindexsing
real*8, allocatable :: elcent(:,:), boundgrain(:,:),             &
    nodex(:,:), nodey(:,:), nodez(:,:),mindist(:),vect(:,:),     &
    sortmindist(:),sortedbgrains(:),uniqu(:),       			 &
    oriensh(:,:),zrot(:,:),slpdirrotate(:,:),   				 &
    grainboundindex(:) ,featureboundsnodes(:,:),                 &
    bxvalue(:),byvalue(:),bzvalue(:),btotal(:,:),bangle(:,:),    &
    normarray(:),normslp(:),minangleval(:),minangleval180(:),    &
    minangleactual(:),rnodes(:,:),vectr(:,:),rdistance(:),       &
    nodesnot(:,:),lxtotal(:,:,:),lxnormin(:,:),lxnorm(:),        &
    vectrnorm(:),lxangle(:,:),minxval0(:,:),minxval180(:,:),     &
    minxval0act(:,:),minxval180act(:,:),minxangle(:),            &
    minxangleact(:),lxnodes(:,:),lxnodesindex(:),                &
    slpdirrotatenorm(:),veca(:,:),vecb(:,:),checkangle(:),       &
    vanorm(:),vbnorm(:),ldistance(:),slpplanrotate(:,:),         &
    slpdirrotatea(:,:),slpplanrotatea(:,:),slpplanrotatenorm(:), &
    slpdirrotateanorm(:),slpplanrotateanorm(:),cosphi(:),cosk(:),&
    lm(:),slpdir(:,:),slpplane(:,:),arraychecker(:)


integer, allocatable :: minindex(:),minang0(:,:),minang180(:,:), &
    minangleactloc(:),boundnodegrainb(:),nodesnotindex(:)
real*8 :: coords(3),vectotal(3),                                 &
    euler1,totalrot(3,3),groundbundid,pi,charctlenght,           &
    twovoxdist(3),centrodes(3),totalrota(3,3)
real :: start, finish,starttotal,finishtotal
 dimension :: SLPDIR1(3,ND), SLPNOR1(3,ND)

LOGICAL, allocatable :: mk(:)
        
logical :: locgrain= .False.



    PI=4.D0*DATAN(1.D0)

 !adjust the slip plane normal and slip direction to removing 
 !any trailing zeros
    allocate(slpdir(nslptl,3))
    allocate(slpplane(nslptl,3))
    
    do i=1,nslptl
        do j=1,3
            slpdir(i,j)=slpdir1(j,i)
            slpplane(i,j)=slpnor1(j,i)
        end do
    end do
 
! check to remove values of interest ** this for loop can be 
! replaced with 'findloc'**

    checkerloop: DO i=1,nodeout
        if(int(nodex(int(feature),i)) == 9999999)then
            arraysize=i-1
			exit checkerloop
		else
		    arraysize=i
			
		end if
     end do checkerloop
     

     
!find the closest boundary.  This will form 'Grain B' from which 
!the information on orientation is extracted to rotate the slip
!systems

!create a vector from the current voxel position to boundary nodes
	allocate (vect(3,arraysize))

    vect(1,:)=nodex(int(feature),1:arraysize)-coords(1)
    vect(2,:)=nodey(int(feature),1:arraysize)-coords(2)
    vect(3,:)=nodez(int(feature),1:arraysize)-coords(3)
    vectotal=(/vect(1,:),vect(2,:),vect(3,:)/)
   	
	allocate (mindist(arraysize))
	allocate (sortmindist(arraysize))
	allocate (mk(arraysize))
	allocate (minindex(arraysize))
	allocate (sortedbgrains(arraysize))
      
! calculate the distance of the current location to the grain boundary and
! multiplying to bring distances to microns
	mindist=(((vect(1,:)**2.D0) + (vect(2,:)**2.D0)+           &
                   (vect(3,:)**2.D0))**0.5D0)
               

!finding the closest boundary grain
	minindexsing=minloc(mindist,1)

!using the index of the minimun distance, the closest grain can be found
!and labelled grain B  
    grainb=int(boundgrain(feature,minindexsing))

    allocate( slpdirrotate(size(slpdir,1),3))
    allocate( slpplanrotate(size(slpplane,1),3))
    
!rotate the slip systems using the angles from each identified grain
  
    call eulercosmatrix(oriensh(grainb,:),totalrot)
    do j =1, size(slpdir,1)
        slpdirrotate(j,:)= matmul(totalrot, slpdir(j,:))
        slpplanrotate(j,:)=matmul(totalrot, slpplane(j,:))
    end do
    
   
    
    allocate( slpdirrotatea(size(slpdir,1),3))
    allocate( slpplanrotatea(size(slpplane,1),3))
    !rotate slip systems for the current grain
    call eulercosmatrix(oriensh(feature,:),totalrota)
    do j =1, size(slpdir,1)
        slpdirrotatea(j,:)= matmul(totalrota, slpdir(j,:))
        slpplanrotatea(j,:)=matmul(totalrota, slpplane(j,:))
    end do
    
  
    
!using the identified closest grain, construct vectors extending from the
!current voxel position to the grain boundary.

!extract only the nodes at the shared boundary between grain A and 
!grain B.


    
!firtly find the size of the array at the indentified grain B
! check to remove values of interest **this for loop can be replaced with 
! 'findloc' if it is available**

    checkerloopb: DO i=1,nodeout
        if(int(nodex(grainb,i)) == 9999999)then
            arraysizeb=i-1
			exit checkerloopb
		else
		    arraysizeb=i
			
		end if
     end do checkerloopb
    
    

    !************************************************************        
    !find index of array locations which are true.
    !For loops are used here to make it more generically applicable
    !but you can also use 'all' first to create a boolean array
    ! and the 'findloc' to find .TRUE. values.
    
  
    
    val=1
    do i=1,arraysizeb
        if(int(boundgrain(grainb,i)) .eq. int(feature))then
            val=val+1
        end if
    end do
    
    allocate (grainboundindex(val-1))
    val=1
    do i=1,arraysizeb
        if(int(boundgrain(grainb,i)) .eq. int(feature))then
            grainboundindex(val)=i
            val=val+1
        end if
    end do
    
  
   !**********************************************************
  

   !find the coordinates of the identified shared boundary nodes
   allocate (featureboundsnodes(val-1,3))
   featureboundsnodes=reshape((/nodex(grainb,grainboundindex), &
         nodey(grainb,grainboundindex),                       &
         nodez(grainb,grainboundindex)/),(/val-1,3/))
    
   !calculate the vector from the boundary nodes to the current
   !voxel
   allocate (bxvalue(val-1))
   allocate (byvalue(val-1))
   allocate (bzvalue(val-1))
   bxvalue=featureboundsnodes(:,1)-coords(1)
   byvalue=featureboundsnodes(:,2)-coords(2)
   bzvalue=featureboundsnodes(:,3)-coords(3)
    
   allocate (btotal(val-1,3)) 
   btotal=reshape((/bxvalue,byvalue,bzvalue/),(/val-1,3/))  
   
   !find the closest angle between the vectors from the interface
   !boundary to the voxel and the slip direction in grain B.
   
   allocate (bangle(size(slpdirrotate,1),val-1))
   
   !these subroutine calls can just be replaced with the 
   !intrisic function norm2*********
   allocate (normarray(val-1))
   allocate (normslp(size(slpdirrotate,1)))
   

   call enorm(btotal,val-1,normarray)
   call enorm(slpdirrotate,size(slpdirrotate,1),normslp)
   !****************************************************
   !determine if the current voxel is on the boundary.
   !If this is true, this voxel should form the grain boundary voxel
   !To do this, determine if any euclidean distances are less than 
   !the elements characteristic length
   
   !calculate characteristic length of voxel
   
   twovoxdist=(elcent(3,:)-elcent(2,:))/2.D0
   charcterlengt=(twovoxdist(1)**2.D0 + twovoxdist(2)**2.D0 +  &
        twovoxdist(3)**2.D0)**0.5D0
   
  !****removed this to see if an improvement appears with the 
  !****voxel at the boundary      
  !this can be done with findloc if available
  do i=1,size(normarray,1)
     if(normarray(i) .lt. charcterlengt)then
        locgrain=.True.
     end if
     
  end do
  !****************** 
   
   allocate (minang0(size(slpdirrotate,1),1))
   allocate (minangleval(size(slpdirrotate,1)))
   allocate (minang180(size(slpdirrotate,1),1))
   allocate (minangleval180(size(slpdirrotate,1)))
   allocate (minangleactual(size(slpdirrotate,1)))
   allocate (minangleactloc(size(slpdirrotate,1)))
   allocate (boundnodegrainb(size(slpdirrotate,1)))
   allocate (rnodes(size(slpdirrotate,1),3))
   allocate (rdistance(size(slpdirrotate,1)))
   
 
   
   if (locgrain == .FALSE.) then
   !voxel not on the boundary


       do i=1,size(slpdirrotate,1)
       !calculate angle between the slipdirection in grian B and the vector
       !created between the current voxel and the boundary nodes.
            bangle(i,:)=dacos(matmul(btotal,slpdirrotate(i,:))/    &
                (normarray*normslp(i))) 
       !find minimum angle and the location index
       !of the array for each slip direction 
            minang0(i,:)=minloc(bangle(i,:))
            minangleval(i)=bangle(i,minang0(i,1))
            
        !find the minimums at angles of 180 degrees
            minang180(i,:)=minloc(abs(bangle(i,:)-pi))
            minangleval180(i)=abs(bangle(i,minang180(i,1))-pi)
        !find the smaller of the 0 and 180 degrees to determine
        !the true mininum value
            if (minangleval(i) .lt. minangleval180(i)) then
                minangleactual(i)=minangleval(i)
                minangleactloc(i)=minang0(i,1)
            else
                minangleactual(i)=minangleval180(i)
                minangleactloc(i)=minang180(i,1)
            end if        
       end do 

       !using the minimum angles, the nodes at the boundary  which connect to the voxel location
       !can be determined
       
       
       
       boundnodegrainb=grainboundindex(minangleactloc)
       rnodes=reshape((/nodex(grainb,int(boundnodegrainb)),        &
                nodey(grainb,int(boundnodegrainb)),                &
                nodez(grainb,int(boundnodegrainb))/),              &
                (/size(slpdirrotate,1),3/)) 
                
       ! the corresponding vector for each slip direction
       allocate (vectr(size(slpdirrotate,1),3))
      
       
       vectr=reshape((/nodex(grainb,int(boundnodegrainb))-coords(1),&
                       nodey(grainb,int(boundnodegrainb))-coords(2),&
                       nodez(grainb,int(boundnodegrainb))-coords(3) &
                       /),(/size(slpdirrotate,1),3/))
                       
       ! calculate the distance from the voxel centroid to the boundary nodes
       call enorm(vectr,size(slpdirrotate,1),rdistance)
   else
   !since the voxel is on the boundary, the distance to the boundary is taken 
   !as half the voxel characteristic length
       
       rdistance(1:size(slpdirrotate,1))=charcterlengt*0.5D0

       do i=1,3
          rnodes(1:12,i)=elcent(int(noel),i)
       end do
   end if
     
   ! calculate the distance from the boundary to boundary in grain B along
   ! the slip direction
   
   !****************************************************************
   !This is to be done by excluding the nodes at the boundary.
   !For loops are used here to make it more generically applicable
   !but you can also use 'all' first to create a boolean array
   ! and the 'findloc' to find .TRUE. values.
   
    valb=1
    do i=1,arraysizeb
        if(int(boundgrain(grainb,i)) .ne. int(feature))then
            valb=valb+1
        end if
    end do
    
    allocate (nodesnotindex(valb-1))
    valb=1
    do i=1,arraysizeb
        if(int(boundgrain(grainb,i)) .ne. int(feature))then
            nodesnotindex(valb)=i
            valb=valb+1
        end if
    end do

    !**********************************************************
    !node coordinates not at the shared boundary.
    allocate (nodesnot(valb-1,3))
    
    nodesnot=reshape((/nodex(grainb,int(nodesnotindex)),      &
                   nodey(grainb,int(nodesnotindex)),          &
                   nodez(grainb,int(nodesnotindex))           &
                   /),(/valb-1,3/))
                   
    
    allocate (lxtotal(size(slpdirrotate,1),valb-1,3))
    allocate (lxnorm(valb-1))
    allocate (lxnormin(valb-1,3))
    allocate (lxangle(size(slpdirrotate,1),valb-1))
    allocate (minxval0(size(slpdirrotate,1),1))
    allocate (minxval180(size(slpdirrotate,1),1))
    allocate (minxval180act(size(slpdirrotate,1),1))
    allocate (minxval0act(size(slpdirrotate,1),1))
    allocate (minxangle(size(slpdirrotate,1)))
    allocate (minxangleact(size(slpdirrotate,1)))
   
    !calculate vector from the boundary to the nodes in Grain B not 
    !at the boundary
    
    ! check to see if the current voxel is at the boundary or not
    ! If not, progress in the following
    allocate (vectrnorm(size(slpdirrotate,1)))
    allocate (slpdirrotatenorm(size(slpdirrotate,1)))
    if (locgrain == .FALSE.)then
     
        call enorm(vectr,size(slpdirrotate,1),vectrnorm)
        do i=1,size(slpdirrotate,1)
            lxtotal(i,:,:)=reshape((/nodex(grainb,int(nodesnotindex))   &
                           -rnodes(i,1),                                &
                           nodey(grainb,int(nodesnotindex))-rnodes(i,2),&
                           nodez(grainb,int(nodesnotindex))-rnodes(i,3) &
                           /),(/valb-1,3/))
            lxnormin=lxtotal(i,:,:)
            !******can be replaced with norm2 function
            call enorm(lxnormin,valb-1,lxnorm)
            !find the angle between the vectors from the shared boundary to 
            !boundary nodes in grain b and the vector developed from the voxel
            !to the shared boundary
            lxangle(i,:)= dacos(matmul(lxtotal(i,:,:),vectr(i,:))/(lxnorm*  &
                        vectrnorm(i)))
            minxval0(i,:)=minloc(lxangle(i,:))
            minxval0act(i,:)=lxangle(i,int(minxval0(i,:)))
            minxval180(i,:)=minloc(abs(lxangle(i,:)-pi))
            minxval180act(i,:)=abs(pi-lxangle(i,int(minxval180(i,:))))
            
            minxangle(i)=minxval0(i,1)
            minxangleact(i)=minxval0act(i,1)
            
         end do
         
     
     else
        
        call enorm(slpdirrotate,size(slpdirrotate,1),slpdirrotatenorm)
        
        do i=1,size(slpdirrotate,1)
            lxtotal(i,:,:)=reshape((/nodex(grainb,int(nodesnotindex))   &
                           -rnodes(i,1),                                &
                           nodey(grainb,int(nodesnotindex))-rnodes(i,2),&
                           nodez(grainb,int(nodesnotindex))-rnodes(i,3) &
                           /),(/valb-1,3/))
            lxnormin=lxtotal(i,:,:)
            !******can be replaced with norm2 function
            call enorm(lxnormin,valb-1,lxnorm)
            !find the angle between the vectors from the shared boundary to 
            !boundary nodes in grain b and the vector developed from the voxel
            !to the shared boundary
            lxangle(i,:)= dacos(matmul(lxtotal(i,:,:),slpdirrotate(i,:))/(lxnorm*  &
                        slpdirrotatenorm(i)))
            minxval0(i,:)=minloc(lxangle(i,:))
            minxval0act(i,:)=lxangle(i,int(minxval0(i,:)))
            minxval180(i,:)=minloc(abs(lxangle(i,:)-pi))
            minxval180act(i,:)=abs(pi-lxangle(i,int(minxval180(i,:))))
            
            !since the vectors in this case can be parallel but in opposite directions
            !an 'if' loop is used all sort through the 0deg and 180deg possible 
            !minimums
            if(minxval0act(i,1)<minxval180act(i,1))then   
                minxangle(i)=minxval0(i,1)
                minxangleact(i)=minxval0act(i,1)
            else
                minxangle(i)=minxval180(i,1)
                minxangleact(i)=minxval180act(i,1)
            end if
            
         end do
     end if
     
     
     allocate (lxnodesindex(size(slpdirrotate,1)))
     allocate (lxnodes(size(slpdirrotate,1),3))
     !using the location of the closest vector in grain B to the developed vector
     !along grain B's slip systems in grain A, the vector and corresponding nodes
     !can be determined.
     lxnodesindex=nodesnotindex(minxangle)
     lxnodes=reshape((/nodex(grainb,lxnodesindex),               &
         nodey(grainb,lxnodesindex),                             &
         nodez(grainb,lxnodesindex)/),(/size(slpdirrotate,1),3/))

     
     !need to check the angle formed by the lxnodes vector and the slip vector
     !in grain A.  If the difference is large (greater than a specified minimum
     !allowable angle), the slip distance in grian B should be ignored since it
     !is unreliable
     allocate (vecb(size(slpdirrotate,1),3))
     allocate (veca(size(slpdirrotate,1),3))
     if (locgrain == .FALSE.)then
        !vector from shared bounday nodes to 
        vecb=lxnodes-rnodes
        !vector from voxel to shared boundary
        do i=1,size(slpdirrotate,1)
            veca(i,:)=elcent(int(noel),:)-rnodes(i,:)
        end do
      else
        !vector from shared bounday nodes to 
        vecb=lxnodes-rnodes
        !voxel at the bounary so use slip direction in grain B
        veca=slpdirrotate
      end if
      
      !calculate the angle between vector a and b
      allocate (vanorm(size(slpdirrotate,1)))
      allocate (vbnorm(size(slpdirrotate,1)))
      allocate (checkangle(size(slpdirrotate,1)))
      allocate (ldistance(size(slpdirrotate,1)))
      
      call enorm(veca,size(slpdirrotate,1),vanorm)
      call enorm(vecb,size(slpdirrotate,1),vbnorm)
      
      !slip distance
      ldistance=vbnorm
      
      
      !calculate the Luster-Morris parameter
      allocate (slpplanrotatenorm(size(slpplanrotate,1)))
      allocate (slpdirrotateanorm(size(slpdir,1)))
      allocate (slpplanrotateanorm(size(slpplane,1)))
      
      call enorm(slpplanrotate,size(slpplanrotate,1),slpplanrotatenorm)
      call enorm(slpdirrotatea,size(slpdir,1),slpdirrotateanorm)
      call enorm(slpplanrotatea,size(slpplane,1),slpplanrotateanorm)
      
      allocate (cosphi(size(slpplane,1)))
      allocate (cosk(size(slpplane,1)))
      allocate (lm(size(slpplane,1)))
     
      do i=1,size(slpplane,1)
        cosphi(i)=dot_product(slpplanrotate(i,:),slpplanrotatea(i,:))/   &
                  (slpplanrotatenorm(i)*slpplanrotateanorm(i))
       
        cosk(i)=dot_product(slpdirrotatea(i,:),slpdirrotate(i,:))/       &
                  (slpdirrotateanorm(i)*normslp(i))
        lm(i)=abs(cosphi(i))*abs(cosk(i))
      end do
      
        

   
     deallocate (vect,mindist  ,sortmindist, mk, minindex,sortedbgrains,  &
                slpdirrotate,slpplanrotate,slpdirrotatea,   &
                slpplanrotatea,grainboundindex,featureboundsnodes,        &
                bxvalue,byvalue,bzvalue,btotal,bangle,normarray,normslp,  &
                minangleval,minang180,minangleval180,minang0,             &
                minangleactual,minangleactloc,boundnodegrainb,rnodes,     &
                nodesnotindex,nodesnot,lxtotal,lxnorm,lxangle,lxnormin,   &
                minxval0,minxval180,minxval180act,minxval0act,minxangle,  &
                minxangleact,vectrnorm,slpdirrotatenorm,lxnodesindex,     &
                lxnodes,vecb,veca,vanorm,vbnorm,checkangle,     &
                slpplanrotatenorm,slpdirrotateanorm,slpplanrotateanorm,   &
                cosphi,cosk)




end subroutine grainsize
!----------------------------------------------------------------------------
! form rotation matrix based on Euler angles
subroutine eulercosmatrix(orient,totalrot)
    Real*8 :: euler1,euler2,euler3,pi,zrot(3,3),xrot(3,3),zrot2(3,3), &
        orient(1,3)
    real*8 :: totalrot(3,3)
    integer :: i, l,j

    PI=4.D0*DATAN(1.D0)
    euler1=orient(1,1)*PI/180.D0
    euler2=orient(1,2)*PI/180.D0
    euler3=orient(1,3)*PI/180.D0
 

    zrot=reshape((/cos(euler1),-sin(euler1),0.d0, &
            sin(euler1),cos(euler1),0.d0, &
                0.d0,0.d0,1.d0/),(/3,3/))
                
    xrot=reshape((/1.d0,0.d0,0.d0, &
            0.d0,cos(euler2),-sin(euler2), &
            0.d0, sin(euler2),cos(euler2)/),(/3,3/))

    zrot2=reshape((/cos(euler3),-sin(euler3),0.d0, &
            sin(euler3),cos(euler3),0.d0, &
                    0.d0,0.d0,1.d0/),(/3,3/))

    !total rotation matrix
    totalrot=transpose(matmul(matmul(zrot2,xrot),zrot))

    end subroutine eulercosmatrix
!------------------------------------------------------------------- 

	!just added this function since norm2 isn't always available
	
	subroutine enorm(array,val,normarray)
	real*8, allocatable, intent(in) :: array(:,:)
	integer*8, intent(in) :: val
	real*8 :: normarray(val)
    

	
	 normarray=(array(:,1)**2.0 + array(:,2)**2.0 + array(:,3)**2.0)**0.50
	
	
	end subroutine enorm



!----------------------------------------------------------------------