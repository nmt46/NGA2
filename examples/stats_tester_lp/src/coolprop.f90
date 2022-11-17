module coolprop
    interface
      function cprop(output,name1,prop1,name2,prop2,fluidname) bind(C,name='PropsSI')
         use iso_c_binding
         real(C_DOUBLE) :: cprop
         character(KIND=c_char), dimension(*) :: output
         character(KIND=c_char), dimension(*) :: name1
         real(C_DOUBLE), value :: prop1
         character(KIND=c_char), dimension(*) :: name2
         real(C_DOUBLE), value :: prop2
         character(kind=c_char), dimension(*) :: fluidname
      end function cprop
!    end interface
! end module coolprop
      function cprop_triv(output,fluidname) bind(C,name='Props1SI')
         use iso_c_binding
         real(C_DOUBLE) :: cprop_triv
         character(KIND=c_char), dimension(*) :: output
         character(KIND=c_char), dimension(*) :: fluidname
      end function cprop_triv
   end interface
end module coolprop
