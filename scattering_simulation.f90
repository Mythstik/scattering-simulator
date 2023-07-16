program particle_scattering
    implicit none

    ! Modify the parameters for the simulation here
    real, parameter :: pi = acos(-1.0)
    real, parameter :: barrier_height = 40.0 ! in eV
    real, parameter :: barrier_width = 0.5 ! in angstroms
    real, parameter :: e_incident = 100.0 ! in eV
    real, parameter :: mass = 1.0 ! in atomic mass units
    real, parameter :: angle_incident = 45.0 ! in degrees
    integer, parameter :: num_particles = 10000
    integer, parameter :: max_scatterings = 1000
    real, parameter :: step_size = 0.01 ! in angstroms
    real, parameter :: angle_tolerance = 0.3 ! in degrees
    character(len=256), parameter :: output_file = "scatterplot.csv"
    
    integer :: i, j
    integer :: num_scatterings
    real :: x, y, z
    real :: vx, vy, vz
    real :: v, e_kinetic, e_potential, e_total
    real :: angle, angle_incident_rad, angle_scattered_rad
    real :: prob_scattering, rand_num
    logical :: scattering_occurred

    call random_seed()
    
    angle_incident_rad = angle_incident * pi / 180.0

    open(unit=10, file=output_file, status="replace")

    write(10, "('x', ',', 'y', ',', 'z', ',', 'angle', ',', 'num_scatterings')")

    do i = 1, num_particles
        x = 0.0
        y = 0.0
        z = 0.0
        vx = sqrt(2.0 * e_incident / mass) * sin(angle_incident_rad)
        vy = 0.0
        vz = sqrt(2.0 * e_incident / mass) * cos(angle_incident_rad)
        num_scatterings = 0
        scattering_occurred = .false.
        
        do j = 1, max_scatterings
            e_potential = barrier_height * step_function(z, barrier_width)
            v = sqrt(vx**2 + vy**2 + vz**2)
            e_kinetic = 0.5 * mass * v**2
            e_total = e_kinetic + e_potential

            prob_scattering = exp(-2.0 * pi * barrier_width / v)

            call random_number(rand_num)
            if (rand_num < prob_scattering) then
                scattering_occurred = .true.
                angle = acos(vz / v) * 180.0 / pi
                angle_scattered_rad = angle * pi / 180.0
                vx = sqrt(2.0 * e_total / mass) * sin(angle_scattered_rad)
                vy = 0.0
                vz = sqrt(2.0 * e_total / mass) * cos(angle_scattered_rad)
                num_scatterings = num_scatterings + 1
            else
                x = x + vx * step_size
                y = y + vy * step_size
                z = z + vz * step_size
            end if
            
            if (z > barrier_width) exit
        end do
        
        if (scattering_occurred) then
            angle = acos(vz / v) * 180.0 / pi
        else
            angle = angle_incident
        end if

        write(10, "(F8.3, ',', F8.3, ',', F8.3, ',', F8.3, ',', I2)") x, y, z, angle, num_scatterings
    end do

    close(10)

contains

    real function step_function(x, width)
        real, intent(in) :: x, width
        if (x < 0.0) then
            step_function = 0.0
        else if (x > width) then
            step_function = 1.0
        else
            step_function = 0.5
        end if
    end function step_function

end program particle_scattering
