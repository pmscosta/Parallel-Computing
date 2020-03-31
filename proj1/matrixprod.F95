program matrixprod
    implicit none
    integer :: lin, col, opt
    logical :: exist
    character(LEN=10000) :: buffer
    print*,"enter two integers"
    CALL get_command_argument(1, buffer)
    read(buffer, *) opt
    CALL get_command_argument(2, buffer)
    read(buffer, *) lin
    CALL get_command_argument(3, buffer)
    read(buffer, *) col
    print *, opt, lin, col
    inquire(file='fortranLog.txt', exist=exist)
    if (exist) then
        open(1, file='fortranLog.txt', status="old", position="append", action="write")
    else
        open(1, file='fortranLog.txt', status="new", action="write")
    end if
    IF (opt == 0) THEN
        call first_implementation(lin, col)
    ELSE IF (opt == 1) THEN
        call second_implementation(lin, col)
    ELSE IF (opt == 2) THEN
        call matmul_implementation(lin, col)
    END IF
end program matrixprod

subroutine first_implementation(m_ar, m_br) 
    implicit none
    real    ::  T1,T2 
    integer :: m_ar, m_br, i, j, k
    real, dimension(m_ar, m_br) :: pha, phb, phc

    do i = 1, m_ar
        do j = 1, m_br
            pha(i, j) = 1.0
        end do
    end do
    
    do i = 1, m_ar
        do j = 1, m_br
            phb(i, j) = i + 1.0
        end do
    end do

    call cpu_time(T1)
    do i = 1, m_ar
        do j = 1, m_br
            do k = 1, m_ar
                phc(i, j) = phc(i,j) + pha(i, k) * phb(k, j)
            end do
        end do
    end do

    call cpu_time(T2)
    print  '("Time = ,",f8.3," seconds.")', T2-T1
    write(1,'(I6,",",I6,",",(f8.3))') m_ar, m_br, T2-T1   
            
end subroutine first_implementation

subroutine second_implementation(m_ar, m_br) 
    implicit none
    real    ::  T1,T2 
    integer :: m_ar, m_br, i, j, k
    real, dimension(m_ar, m_br) :: pha, phb, phc

    do i = 1, m_ar
        do j = 1, m_br
            pha(i, j) = 1.0
        end do
    end do
    
    do i = 1, m_ar
        do j = 1, m_br
            phb(i, j) = i + 1.0
        end do
    end do


    call cpu_time(T1)
    do i = 1, m_ar
        do k = 1, m_ar
            do j = 1, m_br
                phc(i, j) = phc(i,j) + pha(i, k) * phb(k, j)
            end do
        end do
    end do


    call cpu_time(T2)
    print  '("Time = ,",f8.3," seconds.")', T2-T1
    write(1,'(I6,",",I6,",",(f8.3))') m_ar, m_br, T2-T1  
            
end subroutine second_implementation


subroutine matmul_implementation(m_ar, m_br) 
    implicit none
    real    ::  T1,T2 
    integer :: m_ar, m_br, i, j, k
    real, dimension(m_ar, m_br) :: pha, phb, phc

    do i = 1, m_ar
        do j = 1, m_br
            pha(i, j) = 1.0
        end do
    end do
    
    do i = 1, m_ar
        do j = 1, m_br
            phb(i, j) = i + 1.0
        end do
    end do


    call cpu_time(T1)
    
    phc = matmul(pha, phb)

    call cpu_time(T2)
    print  '("Time = ,",f8.3," seconds.")', T2-T1
    write(1,'(I6,",",I6,",",(f8.3))') m_ar, m_br, T2-T1  
            
end subroutine matmul_implementation
