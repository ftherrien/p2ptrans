module JVC_alg

  implicit none

  public :: jvc

contains
  
  subroutine jvc(ifindmax,n,assigncost,absmaxcost,resolution,res,z) 
     !#################################################################
     !>                                                                
     !> brief jonker volgenant castanon method of solving assignment problem
     !>                                                                
     !> details this subroutine solves the full density linear assignment problem
     !>          according to:                                         
     !>       1) "a shortest augmenting path algorithm for dense and sparse linear
     !>          assignment problems," computing 38, 325-340, 1987 by  
     !>          R. Jonker and A. Volgenant, University of Amsterdam.  
     !>                                                                
     !>       2) https://github.com/yongyanghz/lapjv-algorithm-c       
     !>                                                                
     !>       3) https://www.mathworks.com/matlabcentral/fileexchange/26836-
     !>          lapjv-jonker-volgenant-algorithm-for-linear-assignment-v3-0
     !>                                                                
     !> reemark      fixed tolerance problems                          
     !> param[in]   ifindmax    1-find max, 0-find min                 
     !> param[in]   n           order of matrix assigncost             
     !> param[in]   assigncost  cost matrix                            
     !> param[in]   absmaxcost  absolute value of maximum cost (include empty costs)
     !> param[in]   resolution  two reals, differ less then resolution, supposed to be equal
     !> param[out]  res         result of assignment (index of col per row
     !> param[out]  z           summ of assigned costs                 
     !> n                                                              
     !-----------------------------------------------------------------
     ! list of formal parameters                                       
                                             
     integer, intent(in) :: ifindmax ! 1-find max, 0-find min  
                                        
     integer, intent(in) :: n ! order of matrix assigncost   
                                               
     double precision, intent(in) :: assigncost(n,n)  ! cost square matrix   
                             
     double precision, intent(in) :: absmaxcost ! absolute value of maximum cost (include 
                            
     double precision, intent(in) :: resolution  ! two reals, differ less then resolution, supposed to be equal 
                           
     integer, intent(out) :: res(n)   ! result of assignment (index of col per roe
                                   
     double precision, intent(out) :: z   ! summ of assigned costs          
     !-----------------------------------------------------------------
     integer rowsol(n),colsol(n) 
     double precision u(n),v(n) 
     integer i,imin,numfree,prvnumfree,f,i0,k,freerow,pred(n),free1(n),&
    & j,j1,j2,endofpath,last,low,up,loopcnt,collist(n),matches(n)      
     double precision dmin,h,umin,usubmin,v2,d(n) 
     logical unassignedfound 
     double precision cost(n,n) 
                                                                       
     i=0;imin=0;numfree=0;prvnumfree=0;f=0;i0=0;k=0;freerow=0;pred=0 
     free1=0;j=0;j1=0;j2=0;endofpath=0;last=0;low=0;up=0;loopcnt=0 
     collist=0;matches=0 
     dmin=0.;h=0.;umin=0.;usubmin=0.;v2=0.;d=0. 
     unassignedfound=.false. 
     res=0 
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! maximum finding (*-1) - comment this block if minimum is needed 
     !do i=1,n 
     !  do j=1,n 
     !    if(ifindmax.eq.1)then 
                                      ! find max                       
     !      cost(j,i)=-assigncost(j,i) 
     !    elseif(ifindmax.eq.0)then 
                                     ! find min                        
     !      cost(j,i)=assigncost(j,i) 
     !    else 
     !      res=0 
     !      return 
     !    endif 
     !  enddo 
     !enddo 
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! column reduction                                                
                 ! reverse order gives better results                  
     do j=n,1,-1 
       dmin=cost(j,1) 
       imin=1 
       do i=2,n 
         if(cost(j,i).lt.dmin)then 
           dmin=cost(j,i) 
           imin=i 
         endif 
       enddo 
       v(j)=dmin 
       matches(imin)=matches(imin)+1 
       if(matches(imin).eq.1)then 
         ! init assignment if minimum row assigned for first time      
         rowsol(imin)=j 
         colsol(j)=imin 
       else if(v(j).lt.v(rowsol(imin)))then 
         j1=rowsol(imin) 
         rowsol(imin)=j 
         colsol(j)=imin 
                      !0!                                              
         colsol(j1)=-1 
       else 
                     !0! ! row already assigned, column not assigned   
         colsol(j)=-1 
       endif 
     enddo 
                                                                       
     ! reduction transfer                                              
     numfree=0 
     do i=1,n 
                               ! fill list of unassigned 'free' rows   
       if(matches(i).eq.0)then 
         numfree=numfree+1 
         free1(numfree)=i 
                                    ! transfer reduction from rows that are assigned once
       else if(matches(i).eq.1)then 
         j1=rowsol(i) 
         dmin=absmaxcost 
         do j=1,n 
           if(j.ne.j1)then 
             if((cost(j,i)-v(j)).lt.(dmin+resolution))then 
               dmin=(cost(j,i)-v(j)) 
             endif 
           endif 
         enddo 
         v(j1)=v(j1)-dmin 
       endif 
     enddo 
                                                                       
     ! augmenting row reduction                                        
     do loopcnt=1,2 
       ! scan all free rows                                            
       ! in some cases, a free row may be replaced with another one to be scanned next
       k=1 
       prvnumfree=numfree 
       numfree=0 
       do while(k.le.prvnumfree) 
         i=free1(k) 
         k=k+1 
         ! find minimum and second minimum reduced cost over columns   
         umin=(cost(1,i)-v(1)) 
         j1=1 
                           !big                                        
         usubmin=absmaxcost 
         do j=2,n 
           h=(cost(j,i)-v(j)) 
           if(h.lt.usubmin)then 
             if((h.gt.umin).or.(abs(h-umin).lt.resolution))then 
               usubmin=h 
               j2=j 
             else 
               usubmin=umin 
               umin=h 
               j2=j1 
               j1=j 
             endif 
           endif 
         enddo 
         i0=colsol(j1) 
         if((usubmin-umin).gt.resolution)then 
           ! change the reduction of the minimum column to increase the minimum
           ! reduced cost in the row to the subminimum                 
           v(j1)=v(j1)-(usubmin-umin) 
              ! minimum and subminimum equal                           
         else 
                            ! minimum column j1 is assigned            
           if(i0.gt.-1)then 
             ! swap columns j1 and j2, as j2 may be unassigned         
             j1=j2 
             i0=colsol(j2) 
           endif 
         endif 
         ! (re-)assign i to j1, possibly de-assigning an i0            
         rowsol(i)=j1 
         colsol(j1)=i 
         if(i0.gt.-1)then 
           if((usubmin-umin).gt.resolution)then 
             ! put in current k, and go back to that k                 
             ! continue augmenting path i - j1 with i0                 
             k=k-1 
             free1(k)=i0 
           else 
             ! no further augmenting reduction possible                
             ! store i0 in list of free rows for next phase            
             numfree=numfree+1 
             free1(numfree)=i0 
           endif 
         endif 
       enddo 
     enddo 
                                                                       
     ! augment solution for each free row                              
     do f=1,numfree 
                        ! start row of augmenting path                 
       freerow=free1(f) 
       ! dijkstra shortest path algorithm                              
       ! runs until unassigned column added to shortest path tree      
                !do j=n,1,-1                                           
       do j=1,n 
         d(j)=cost(j,freerow)-v(j) 
         pred(j)=freerow 
                      ! init column list                               
         collist(j)=j 
       enddo 
       ! columns in 0..low-1 are ready, now none.                      
       ! columns in low..up-1 are to be scanned for current minimum, now none
       ! columns in up..dim-1 are to be considered later to find new minimum
       ! at this stage the list simply contains all columns            
       low=1 
       up=1 
       unassignedfound=.false. 
       do while(.not.unassignedfound) 
                           ! no more columns to be scanned for current minimum
         if(up.eq.low)then 
           last=low-1 
           ! scan columns for up..dim-1 to find all indices for which new minimum occurs
           ! store these indices between low..up-1 (increasing up)     
           dmin=d(collist(up)) 
           up=up+1 
           do k=up,n 
             j=collist(k) 
             h=d(j) 
             if((h.lt.dmin).or.(abs(h-dmin).lt.resolution))then 
                                 ! new minimum                         
               if(h.lt.dmin)then 
                 up=low 
                 dmin=h 
               endif 
               ! new index with same minimum, put on undex up, and extend list
               collist(k)=collist(up) 
               collist(up)=j 
               up=up+1 
             endif 
           enddo 
           ! check if any of the minimum columns happens to be unassigned
           ! if so, we have an augmenting path right away              
           do k=low,(up-1) 
                                             !                         
             if(colsol(collist(k)).lt.1)then 
               endofpath=collist(k) 
               unassignedfound=.true. 
                    ! break do loop                                    
               exit 
             endif 
           enddo 
         endif 
         if(.not.unassignedfound)then 
           ! update 'distances' between freerow and all unscanned columns, via next scanned column
           j1=collist(low) 
           low=low+1 
           i=colsol(j1) 
           h=(cost(j1,i)-v(j1)-dmin) 
           do k=up,n 
             j=collist(k) 
             v2=(cost(j,i)-v(j)-h) 
             if(v2.lt.d(j))then 
               pred(j)=i 
                                                  ! new column found at same minimum value!if(abs(v2-dim).lt.eps)then ! new column found at same minimum value
               if(abs(v2-dmin).lt.resolution)then 
                 if(colsol(j).lt.0)then 
                   ! if unassigned, shortest augmenting path is complete
                   endofpath=j 
                   unassignedfound=.true. 
                        ! break do loop                                
                   exit 
                 else 
                   ! else add to list to be scanned right away         
                   collist(k)=collist(up) 
                   collist(up)=j 
                   up=up+1 
                 endif 
               endif 
               d(j)=v2 
             endif 
           enddo 
         endif 
       enddo 
                                                                       
       ! update column prices                                          
       do k=1,last 
         j1=collist(k) 
         v(j1)=v(j1)+d(j1)-dmin 
       enddo 
                                                                       
       ! reset row and column assignments along the alternating path   
       do while(i.ne.freerow) 
         i=pred(endofpath) 
         colsol(endofpath)=i 
         j1=endofpath 
         endofpath=rowsol(i) 
         rowsol(i)=j1 
       enddo 
     enddo 
                                                                       
     ! calculate optimal cost                                          
     z=0.0 
              ! do i=n,1,-1                                            
     do i=1,n 
       j=rowsol(i) 
       u(i)=(cost(j,i)-v(j)) 
       z=z+assigncost(j,i) 
                        ! form result                                  
       res(i)=rowsol(i)!colsol(i) 
     enddo 
     print *, res                                                                  
  end subroutine
      
end module JVC_alg                                       
