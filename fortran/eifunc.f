	  implicit real*8 (a-h, o-z)
           
             write *, 'put x' 
               read *, x
			   call eix(f,ei)
			   s=0.d0
			   do i=1,x
			   s=s+ f
			   end do
			   write *, 'no of interval', x  , 'result is',s
             end
              
		function f(x)
             implicit real*8 (a-h, o-z)
             f=x
             end 
			   

			