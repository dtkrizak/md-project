Module mdproj
!*************************************************************************
! The Module initializes the following variables for the program:
!   fij    :  sum of forces exerted on particle i from particles j = [1,N]
!   tstep  :  time step
!   en     :  total energy of system
!   temp   :  temperature
!   press  :  pressure
!   rdist  :  distance between two particles
!   kb     :  boltzmann constant
!   sigma  :  diameter of particle
!   num    :  number of particles
!   dens   :  density of the particles
!   vol    :  volume of the box the particles are in (assumed its a cube)
!   vel    :  velocity of the particles (vx, vy, vz)
!   temp0  :  inital temperature estimate to calculate v0
!   mass   :  mass of a particle
!   pi     :  define pi (3.14159)
!   posit  :  position of particles (x, y, z)
!   deltap :  distance between particles i and j
!   force  :  force exerted on particle i from particle j
!   fx     :  force in the x-direction
!   fy     :  force in the y-direction
!   fz     :  force in the z-direction
!   cx     :  difference of distance of particle i and j in x-direction
!   cy     :  difference of distance of particle i and j in y-direction
!   cz     :  difference of distance of particle i and j in z-direction
!	rmsvel :  the root mean squared velocity
!
! Input:   dens
! Output:  en, temp, press, vol
!
! for an inert gas such as Ar:
!   mass = 6.69E-26 kg
!   sigma = 0.34 nm = 0.34E-9
!	http://www.sklogwiki.org/SklogWiki/index.php/Argon
!
! for Xe:
!	mass = 2.1801714E-25 kg
!	sigma = 4.1E-10 m 
!	http://chem.libretexts.org/Core/Physical_and_Theoretical_Chemistry/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Intermolecular_Forces/Specific_Interactions/Lennard-Jones_Potential
!
! The units utilized to carry out the calculation are based on the 
! relative length of an atom, the mass of an atom, and kbT
!	mass = 1.
!	kbT = 1.
!	sigma = 1.
!	cutoff = 2.5*sigma
!*************************************************************************

integer, parameter :: DP = SELECTED_REAL_KIND(8)
!real, parameter :: tstep = 1.E-10
!Sets accuracy, time step and max time
real (Kind = DP), parameter:: tmax=0.1, tstep = 1.E-5, accuracy = 1.E-25
real (Kind = DP) :: en, temp, temp2, press, press2, dens, vol
real (Kind = DP) :: rmsvel, rmsvel2, en2, v_frac
!real (Kind = DP) :: temp0, equilT, equilE, equilP, equilVel
real (Kind = SELECTED_REAL_KIND(8)) :: side
real, parameter :: cutoff=2.5, kbT = 2.
!real, parameter :: kb=1.38064852E-23, h = 6.62607004E-34
integer, parameter :: num = 300
!real :: maxcalc
real (Kind = DP), dimension(num, 3) :: vel0, vel, vel2, posit0, posit, posit2
real (Kind = DP), dimension(num, 3) :: vel_NVE, vel_And, posit_NVE, posit_And
real (Kind = DP), dimension(num, num):: deltap, deltap2
real (Kind = DP), dimension(num, 3) :: fij, fij2
real (Kind = DP), dimension(3):: rms_dis_NVE, rms_dis_And
real, parameter :: pi = 4. * atan(1.), sigma = 1., v_frac_unit = 1.0E-1
real, parameter :: mass = 1.
real (Kind = DP):: t, velavg, velmin, velmax, velstddev
real, parameter :: half = 0.5
!real, parameter :: x = 1E10
!integer :: acc_count
!v_frac is # change per timestep
integer:: counter1

End Module mdproj


Program project
!*************************************************************************
! This program performs a NVE molecular dynamics simulation using
! the velocity verlet algorithm
!*************************************************************************
Use mdproj

dens = 0.8
!Read(*,*), dens
vol = num / dens
side = vol**(1./3.)

!Print *, 'At this density, with N =', num, ' particles, the volume is', vol
Print *, 'side', side

!Next the starting velocity of the particles are estimated using
! the Boltzmann distribution given a starting temperature
!EDIT: Currently, the velocities are calculated using kbT, not temp0

!Print *, 'Enter an estimate of the temperature (K):'
!temp0 = 100
!Read(*,*), temp0
!Print *, 'The estimated temperature is', temp0, 'K'

!Initializes the velocities
call init_posit(posit0)
call init_v0(vel0)
!Prints initial position and velocity
!Print *, ' '

!Do i = 1, num
!	posit0()
!Print *, 'posit0', posit0
!Print *, ' '
!Print *, 'vel0', vel0


!Sets the accuracy of the calculations for setting values equal to zero
!acc_count=0

call vel_rms(vel, rmsvel)
!Print *, ' '
!Print *, 'rmsvel', rmsvel

!Calculate postion (posit) and velocity (vel) using velocity verlet
! approximation which continues to iterate until time, t, has 
! reached max time, tmax.
vel = vel0
posit = posit0
call NVE_sim()

!Stores the values of velocity and position for the NVE simulation
vel_NVE = vel2
posit_NVE = posit2

!Prints final values of NVE sim
!Print *, 'side:', side
Print *, ' '
!Print *, 'posit_NVE', posit_NVE
!Print *, ' '
Print *, 'vel_NVE', vel_NVE
!Print *, ' '
Print *, 'temp2', temp2
Print *, 'press2', press2
Print *, 'en2', en2
Print *, 'rmsvel2', rmsvel2

!Calculates the root mean square dispalcement for the NVE Simulation
Call displace_rms(rms_dis_NVE)
Print *, 'NVE rms_dis', rms_dis_NVE

!Calculates the average, minimum, maximum, and standard deviation for the velocity
Call vel_analysis(vel_NVE, velavg, velmin, velmax, velstddev)
Print *, 'NVE: velavg', velavg, 'velmin', velmin, 'velmax', velmax, 'velstddev', velstddev

!Calculates the rms displacement for the simulation at the final time step
! compared to the initial position
!Call displace_rms(rms_dis_NVE)
!Print *, 'NVE rms_dis', rms_dis_NVE

!Calculate postion (posit) and velocity (vel) using velocity verlet
! approximation with an Andersen thermostat, which continues to
! iterate until time, t, has reached max time, tmax.
!
!The velocity array, vel, is set to the initial condition calculated
! previously. The position array, posit, is recalculated based on
! the size of the box
!
!The density is changed and to determine the different pressure
! at densities between 0.1 and 0.8 with increments of 0.1
!
!v_frac is # of part that change per timestep for the Andersen thermostat
v_frac = v_frac_unit * tstep

Do counter1 = 1, 8
!Do counter1 = 1, 2
	dens = Real(counter1, 8) * 0.1
!	dens = 0.8
!	Print *, 'counter1', counter1, 'density', dens
	vol = num / dens
	side = vol**(1./3.)
	vel = vel0
	call init_posit(posit)
	call Andersen_sim()
	Print *, 'dens =', dens,'  side:', side
	!Stores the values of velocity and position for the Andersen simulation
	vel_And = vel2
	posit_And = posit2

!	Print *, 'side:', side
!	Print *, ' '
!	Print *, 'posit_And', posit_And
!	Print *, ' '
	Print *, 'vel_And', vel_And
!	Print *, ' '
	Print *, 'temp2', temp2
	Print *, 'press2', press2
	Print *, 'en2', en2
	Print *, 'rmsvel2', rmsvel2

	!Calculates the rms displacement for the simulation at the final time step
	! compared to the initial position
	Call displace_rms(rms_dis_And)
	Print *, 'And rms_dis', rms_dis_And

	!Calculates the average, minimum, maximum, and standard deviation for the velocity
	Call vel_analysis(vel_And, velavg, velmin, velmax, velstddev)
	Print *, 'And: velavg', velavg, 'velmin', velmin, 'velmax', velmax, 'velstddev', velstddev

End Do

End Program project



!The Subroutine init_posit initializes the initial positions of the
! particles into a lattice like structure
Subroutine init_posit(positin)
Use mdproj
real :: lattice, latticepart, posit_max, posit_min
, dimension(num, 3)::positin
integer :: dovarx, dovary, dovarz, partcount, latticepartint

posit_max=-1.0
posit_min=100000.0


latticepart = num**(0.33333333333333)
latticepartint = nint(latticepart)
lattice = side/(latticepart)
!Print *, 'latticepart', latticepart
!Print *, 'latticepartint', latticepartint
!Print *, 'lattice', lattice
!STOP
partcount = 1

!Do loop populates posit0
Do dovarz = 0, (latticepartint-1)
	Do dovary = 0, (latticepartint-1)
		Do dovarx = 0, (latticepartint-1)
		if(partcount > num) then
			Exit
		else
			positin(partcount, 1) = lattice * dovarx
			positin(partcount, 2) = lattice * dovary
			positin(partcount, 3) = lattice * dovarz
!			Print *, 'posit0', positin(partcount,1), positin(partcount,2), positin(partcount,3)
!			Print *, 'dovarx, y, z', dovarx, dovary, dovarz
!			Print *, 'partcount', partcount
			partcount = partcount + 1
		End if
		End Do
	End Do
End Do
!Print *, 'dovarx, y, z', dovarx, dovary, dovarz
!Print *, 'partcount', partcount
End Subroutine init_posit


!The subroutine init_v0 sets the inital velocities of the particles
! utilizing the Boltzmann distribution by multiplying the random
! Gaussian generator by sqrt(kbT/mass)
!Reference: https://en.wikibooks.org/wiki/Fortran/Fortran_procedures_and_functions
Subroutine init_v0(velin)
Use mdproj
real :: mean, stdDev, u, v, s, spare, mul, outvar
integer :: dovar1, dovar2
real (Kind = DP), dimension(num, 3):: velin
logical :: isSpareReady = .FALSE.

mean = 0.
stdDev = 1.

!Do loops populate the entire array for vel0
Do dovar1 = 1, num
    Do dovar2 = 1, 3
    isSpareReady = .FALSE.

    !If statement handles generation of the boltmann distribution
    ! for the initial veloities
    If (isSpareReady) then
    isSpareReady = .FALSE.
    outvar = spare * stdDev + mean
    Else
    Do
        Call random_number(u)
        u = u * 2 - 1
        Call random_number(v)
        v = v * 2 - 1
        s = u*u + v*v
        If((s < 1.) .AND. (s /= 0)) then
        Exit
        End if
    End Do
    mul = sqrt(-2.0*log(s)/s)
    spare = v * mul
    isSpareREady=.TRUE.
    outvar = mean + stdDev * u * mul
    End if

	velin(dovar1,dovar2)=sqrt(kbT/mass)*outvar

	End Do
End Do

End Subroutine init_v0


!Loop Calculated the force between the particles by first calculating the distance
! by taking the difference in the x, y, and z directions, calculating the length
! of the vector and storing it in array deltap, then normalizing the vector c to use it
! to break down the force vector into fx, fy, fz and summing the forces exerted on
! particle i in the x(1), y(2), and z(3) directions in array fij
Subroutine force_calc(positin, fij_out, deltap1)
Use mdproj
real (KIND = DP), dimension(num, 3) :: positin, fij_out
real (Kind = DP), dimension(num, num) :: deltap1
real (Kind = DP):: force, fx, fy, fz, cx, cy, cz
integer :: count1, count2, i, j
real :: cmin, cmax
!logical :: truthtest
!infinite = HUGE(infinite)

Do count1 = 1, num
	Do count2 = 1, 3
		fij_out(count1,count2)=0.
	End Do
End Do

cmin = 10
cmax = -10

!Print *, '0positin', positin

Do i = 1, num
	Do j = 1, num
		if(i == j) then
!			Print *, "Entered the if at i=", i, 'and j=', j
			deltap1(i,j) = 0.
!			Print *, '1positin', positin(i,1), 'j', positin(j,1)
		else
!			Print *, '1posit', posit, 'posit2', posit2
!			Print *, '2positin', positin(i,1), 'j', positin(j,1)
!			Print *, "Entered the else if at i=", i, 'and j=', j

			cx = positin(i,1) - positin(j,1)
			cy = positin(i,2) - positin(j,2)
			cz = positin(i,3) - positin(j,3)
!			Print *, 'cx', cx, 'cy', cy, 'cz', cz

!			if(cx > (side/2)) then
!				cx = cx - (side/2)
!			Print *,'cx', cx, 'i', i, 'j', j
!			End if
!			if(cx < (-side/2)) then
!				cx = cx + (side/2)
!			Print *,'cx', cx, 'i', i, 'j', j
!			End if
!
!			if(cy > (side/2)) then
!				cy = cy - (side/2)
!				Print *, 'cy', cy, 'i', i, 'j', j
!			End if
!			if(cy < (-side/2)) then
!				cy = cy + (side/2)
!!				Print *, 'cy', cy, 'i', i, 'j', j
!			End if

!			if(cz > (side/2)) then
!				cz = cz - (side/2)
!				Print *, 'cz', cz, 'i', i, 'j', j
!			End if
!			if(cz < (-side/2)) then
!				cz = cz + (side/2)
!				Print *, 'cz', cz, 'i', i, 'j', j
!			End if

			if(cx < cmin) cmin = cx
			if(cy < cmin) cmin = cy
			if(cz < cmin) cmin = cz
			if(cx > cmax) cmax = cx
			if(cy > cmax) cmax = cy
			if(cz > cmax) cmax = cz

!			Print *, '3positin', positin(i,1), 'j', positin(j,1), 'i,j', i, j
			deltap1(i,j) = sqrt((cx)**2 + (cy)**2 + (cz)**2)
!			Print *, 'deltap1', deltap1(i,j), 'i', i, 'j', j
!			truthtest = (deltap1(i,j) == 0.) .OR. (deltap1(i,j) >= cutoff)
!			Print *, 'Truthtest', truthtest

!				if(deltap(i,j) /= deltap(i,j)) then
!				Print *, 'FATAL ERROR (NaN) at ', i, j
!				End if

!				if(deltap(i,j) >= infinite) then
!				Print *, ' FATAL ERROR (Infinity) at', i, j
!				End if

			if((deltap1(i,j) == 0.) .OR. (deltap1(i,j) >= cutoff)) then
!			if(deltap1(i,j) >= cutoff) then
				fij_out(i,1) = fij_out(i,1) + 0.0
				fij_out(i,2) = fij_out(i,2) + 0.0
				fij_out(i,3) = fij_out(i,3) + 0.0
!				Print *, "the program exited the else if at i=", i, 'and j=', j
!				Print *, 'deltap1', deltap1(i,j), 'i', i, 'j', j
!				Print *, 'positin: i', positin(i,1), 'j', positin(j,1)
!				Print *, 'positin: i', positin(i,2), 'j', positin(j,2)
!				Print *, 'positin: i', positin(i,3), 'j', positin(j,3)
!				Print *, 'cx', cx, 'cy', cy, 'cz', cz
			else
!				Print *, "the program entered the else at i=", i, 'and j=', j
				cx = cx/deltap1(i,j)
				cy = cy/deltap1(i,j)
				cz = cz/deltap1(i,j)

!				if(cx < accuracy) then
!					cx = 0.
!				End if
!				if(cy < accuracy) then
!					cy = 0.
!				End if
!				if(cz < accuracy) then
!					cz = 0.
!				End if

!				Print *, 'cx, y, z', cx, cy, cz, deltap(i,j)
!				force = -1/kbT*((sigma/deltap(i,j))**12 - (sigma/deltap(i,j))**6)/deltap(i,j)
!				Print *, ' force:', force
				force = 4/kbT*((12*sigma**12/(deltap1(i,j)**13)) - (6*sigma**6/(deltap1(i,j)**7)))
				fx = force * cx
				fy = force * cy
				fz = force * cz
				fij_out(i,1) = fij_out(i,1) + fx
				fij_out(i,2) = fij_out(i,2) + fy
				fij_out(i,3) = fij_out(i,3) + fz
!				Print *, 'force', force, 'dist', deltap1(i,j)
!				Print *, 'i', i, 'j', j
!				Print *, 'fx', fx, 'fy', fy, 'fz', fz
!				Print *, '     The sum forces for particle i are:', fij_out(i,1), fij_out(i,2), fij_out(i,3)
!				Print *, ' '
!				Print *, ' deltap', deltap(i,j)
			End if
!				if(fx /= fx) then
!				Print *, 'FATAL ERROR (NaN) fx at ', i, j
!				End if

!				if(fx >= infinite) then
!				Print *, ' FATAL ERROR (Infinity) fx at', i, j
!				End if

!				if(fy /= fy) then
!				Print *, 'FATAL ERROR (NaN) fy at ', i, j
!				End if

!				if(fy >= infinite) then
!				Print *, ' FATAL ERROR (Infinity) fy at', i, j
!				End if

!				if(fz /= fz) then
!				Print *, 'FATAL ERROR (NaN) fz at ', i, j
!				End if

!				if(fz >= infinite) then
!				Print *, ' FATAL ERROR (Infinity) fz at', i, j
!				End if

		End if
!		Print *, "the program completed the if statement at i=", i, 'and j=', j
	End Do
!	Print *, 'i', i, 'fij_outx', fij_out(i,1), 'fij_outy', fij_out(i,2), 'fij_outz', fij_out(i,3)

!	if(fij_outin(i,1) < accuracy) then
!	fij_outin(i,1) = 0.
!	acc_count = acc_count + 1
!	End if
!	If(fij_outin(i,2) < accuracy) then
!	fij_outin(i,2) = 0.
!	acc_count = acc_count + 1
!	End if
!	if(fij_outin(i,3) < accuracy) then
!	fij_outin(i,3) = 0.
!	acc_count = acc_count + 1
!	End if

End Do

!Print *, 'cmin', cmin, 'cmax', cmax
!Print *, '2posit', posit, 'posit2', posit2

End Subroutine force_calc


!The NVE_sim subroutine performs the verlet approximation on
! the position (posit) and velocity (vel) vectors
Subroutine NVE_sim()
Use mdproj
integer :: m, j1, j2, j3, l1
REAL(Kind = DP):: posit_max, posit_min
!logical :: equilTtruth = .FALSE.
!real (Kind = SELECTED_REAL_KIND(16)):: j5
!integer (Kind = Selected_Int_Kind(32))::j4
!infinite = HUGE(infinite)

!Print *, '1posit', posit, 'posit2', posit2
!Sets counter variables
m = 1
t = 0.

posit_max=-1
posit_min=1000000

Do
!	Print *, 'm NVE', m
	!Calculate distance between particles and the forces exerted on the particles
	Call force_calc(posit, fij, deltap)
!	Print *, 'forcecalc ', m, ' completed'
!	Print *, '2posit', posit, 'posit2', posit2
!	posit2 = posit + vel*tstep + (fij/2/mass*(tstep**2))
!	Print *, '3posit', posit, 'posit2', posit2
!	Print *, 'acc_count', acc_count

	Do j1 = 1, num
		Do j2 = 1, 3
!			if(posit(j1,j2) > 10*side) then
!			Print *, 'Error - particle1'
!			Exit
!			Endif

!			if(posit2(j1,j2) > 10*side) then
!			Print *, 'Error - particle2'
!			Exit
!			Endif
			posit2(j1,j2) = posit(j1,j2) + vel(j1,j2)*tstep + (fij(j1,j2)/2/mass*(tstep**2))
!			Print *, '   posit2', posit2(j1,j2), 'posit', posit(j1,j2), 'vel', vel(j1,j2), 'fij', fij(j1,j2) !'tstep', tstep

!			if((vel(j1,j2)*tstep) /= (vel(j1,j2)*tstep)) then
!			Print *, 'FATAL ERROR (NaN) for vel term', vel(j1, j2), tstep
!			End if

!			if((vel(j1,j2)*tstep) >= infinite) then
!			Print *, ' FATAL ERROR (Infinity) for vel term', vel(j1, j2), tstep
!			End if

!			if((fij(j1,j2)/2/mass*(tstep**2)) /= (fij(j1,j2)/2/mass*(tstep**2)) ) then
!			Print *, 'FATAL ERROR (NaN) for fij term', fij(j1, j2), 1/2/mass*(tstep**2)
!			End if

!			if((fij(j1,j2)/2/mass*(tstep**2)) >= infinite) then
!			Print *, ' FATAL ERROR (Infinity) fij term', fij(j1, j2), 1/2/mass*(tstep**2)
!			End if

			if(posit2(j1,j2) >= side) then
				j3 = posit2(j1,j2)/side
				posit2(j1,j2) = posit2(j1,j2) - (j3 * side)
!				if(j3>2) then
!				Print *, 'j3', j3
!				End if
!				Print *, 'j1', j1, 'j2', j2, 'posit2', posit2(j1,j2)
				!Based on the bulk system, the velocity does not change direction
!				vel(j1,j2) = - vel(j1,j2)
!				Print *, 'neg vel', vel(j1, j2)
!				elseif(posit2(j1,j2) == side) then
!				vel(j1,j2) = -vel(j1, j2)
!				Print *, 'negvel', vel(j1, j2)
			End if
			if (posit2(j1,j2) < 0.0) then
				j3 = -(posit2(j1,j2)/side) + 1
				posit2(j1,j2) = posit2(j1,j2) + (j3 * side)
!				if(j3>2) then
!				Print *, 'j3', j3
!				End if
!				Print *, 'neg fix', posit2(j1,j2), j3, (j3*side)
			End if

			IF (posit2(j1,j2)>posit_max) posit_max=posit2(j1,j2)
			IF (posit2(j1,j2)<posit_min) posit_min=posit2(j1,j2)

!			if (posit2(j1,j2) < -0.0001) then
!				posit2(j1,j2) = abs(posit2(j1,j2))
!				j4 = ((posit2(j1,j2)-side)/side)
!				j5 = j4
!				Print *, 'neg fix', posit2(j1,j2), j5, (side*j5)
!				posit2(j1,j2) = posit2(j1,j2) - (side*j5)
!				Print *, 'neg fix pt 2', posit2(j1,j2)
!				Print *, 'j1', j1, 'j2', j2, 'posit2 neg', posit2(j1,j2), 'j3', j3
!			End if

!			Print *, posit2(j1,j2)

		End Do
	End Do
!PRINT*, posit_max, posit_min
!STOP


	!Calculates the force after the first time step, used to calculate vel2
	Call force_calc(posit2, fij2, deltap2)
!	Print *, 'forcecalc2 ', m, ' completed'
!	Print *, 'acc_count', acc_count

	!Calculates the velocity vel2, after the first time step
	Do l1 = 1, num
		vel2(l1, 1) = vel(l1, 1) + ((fij(l1, 1)+fij2(l1, 1))/2/mass*tstep)
		vel2(l1, 2) = vel(l1, 2) + ((fij(l1, 2)+fij2(l1, 2))/2/mass*tstep)
		vel2(l1, 3) = vel(l1, 3) + ((fij(l1, 3)+fij2(l1, 3))/2/mass*tstep)
!		Print *, '   vel2x', vel2(l1,1), 'velx', vel(l1,1), 'fij2x', fij2(l1,1), 'fijx', fij(l1,1)
!		Print *, '   vel2y', vel2(l1,2), 'vely', vel(l1,2), 'fij2y', fij2(l1,2), 'fijy', fij(l1,2)
!		Print *, '   vel2z', vel2(l1,3), 'velz', vel(l1,3), 'fij2z', fij2(l1,3), 'fijz', fij(l1,3)

!		If(vel2(l1,1) < accuracy) then
!			vel2(l1,1) = 0.
!			Print *, 'min vel2 x'
!			acc_count = acc_count + 1
!		End if
!		If(vel2(l1,2) < accuracy) then
!			vel2(l1,2) = 0.
!			Print *, 'min vel2 y'
!			acc_count = acc_count + 1
!		End if
!		If(vel2(l1,3) < accuracy) then
!			vel2(l1,3) = 0.
!			Print *, 'min vel2 z'
!			acc_count = acc_count + 1
!		End if

!		Print *, 'vel', vel(l1, 1), vel(l1, 2), vel(l1, 3)
!		Print *, 'vel2', vel2(l1, 1), vel2(l1, 2), vel2(l1, 3)
!		Print *, 'fij', fij(l1, 1), fij(l1, 2), fij(l1, 3)
!		Print *, 'fij2', fij2(l1, 1), fij2(l1, 2), fij2(l1, 3)
	End Do


	Call vel_rms(vel, rmsvel)
!	Print *, 'The rms vel =', rmsvel
	Call vel_rms(vel2, rmsvel2)
!	Print *, 'The rms vel2 =', rmsvel2

	!Calculates the thermodynamic properties for the initialized variables
	Call thermo_calc(rmsvel, temp, press, en, fij, deltap, posit)
!	Print *, 'rmsvel =', rmsvel
!	Print *, 'temp =', temp
!	Print *, 'press =', press
!	Print *, 'en =', en

	!Calculates the thermodynamic properties for the initialized variables
	Call thermo_calc(rmsvel2, temp2, press2, en2, fij2, deltap2, posit2)
!	Print *, 'rmsvel2 =', rmsvel2
!	Print *, 'temp2 =', temp2
!	Print *, 'press2 =', press2
!	Print *, 'en2 =', en2

	!The following 'equilibrium' variables calculates the difference
	! between the variables: temp, press, en, and rmsvel for before
	! the time step and after the time step
!	equilT = (temp2 - temp)/temp2
!	Print *, '   The equilT =', equilT
!	equilP = (press2 - press)/press2
!	Print *, '   The equilP =', equilP
!	equilE = (en2 - en)/en2
!	Print *, '   The equilE =', equilE
!	equilVel = (rmsvel2 - rmsvel)/rmsvel2
!	Print *, '   The equilVel =', equilVel

	if(t >= tmax) then
		Exit
	else
		!sets vel and posit values for next time step
		vel = vel2
		posit = posit2
		!adds to the counter for the number of iterations completed
		m = m + 1
		t = t + tstep
		Cycle
	End if
End Do


End Subroutine NVE_sim


!Subroutine vel_rms calculates the root mean square velocity of the given velocity
! array (velarray) and outputs the velocity in the variable velout
Subroutine vel_rms(velarray, velout)
Use mdproj
real (Kind = DP), Dimension(num, 3) :: velarray
real (Kind = DP), Dimension(num) :: magarray
real (Kind = DP) :: velout, velsum
integer :: loop1

velsum = 0.

Do loop1 = 1, num
	!Takes the magnitude of the velocity vector of particle i = loop1
	magarray(loop1) = sqrt(velarray(loop1,1)**2 + velarray(loop1,2)**2 + velarray(loop1,3)**2)
!	Print *, 'loop1', loop1
!	Print *, 'velarray', velarray(loop1,1), velarray(loop1,2), velarray(loop1,3)

	!Adds the square of the magnitude of particle i = loop1
	velsum = velsum + magarray(loop1)**2
!	Print *, 'velsum', velsum
End Do

!Print *, 'velsum', velsum, 'magarray', magarray, 'velout', velout
!Outputs the root mean square velocity in the velout variable
velout = sqrt(velsum/num)

End Subroutine vel_rms


!Subroutine thermo_calc takes the root mean square velocity (rms_velin),
! and returns the temperature (tempout), pressure (pressout), 
! and the total energy (enout)
Subroutine thermo_calc(rms_velin, tempout, pressout, enout, fijin, deltapin, positin)
Use mdproj
real (Kind = DP), intent(in):: rms_velin
real (Kind = DP), intent(out)::tempout, pressout, enout
real (Kind = DP) :: ensum, press_force_term
real (Kind = DP), dimension(num, 3), intent(in)::fijin, positin
real (Kind = DP), dimension(num, 3)::Uij
real (Kind = DP), dimension(num, num), intent(in)::deltapin
real (Kind = DP):: force, cx, cy, cz
real (Kind = DP), dimension(num, num, 3)::Ucalc
integer:: e1, d1

!Calculates the temperature based on the root mean squared velocity
tempout = (rms_velin**2) * mass / 3

!Calculates the pressure
pressout = num * tempout / vol
press_force_term = 0.
Do d1 = 1, num
		press_force_term = press_force_term + fijin(d1, 1) + fijin(d1, 2) + fijin(d1, 3)
End Do
press_force_term = press_force_term / vol
pressout = pressout + press_force_term

!Calculates the energy based on ideal gas
Call Uij_calc(positin, Uij, deltapin)
ensum = 0.
Do e1 = 1, num
	ensum = ensum + sqrt(Uij(e1,1)**2 + Uij(e1,2)**2 + Uij(e1,3)**2)
	ensum = ensum + 1/2*mass*sqrt(vel2(e1,1)**2 + vel2(e1,2)**2 + vel2(e1,3)**2)
End Do
enout = ensum
!enout = (rms_velin**2)*mass*num/2 + ensum

End Subroutine thermo_calc



!This subroutine runs a NVT simulation utilizing the Andersen thermostat
! to progress the simulation while holding temperature constant
! by continuously scaling the velocities after each time step
! that is calculated by the velocity position verlet approximation
Subroutine Andersen_sim()
Use mdproj
integer :: m, j1, j2, j3
!logical :: equilPtruth = .FALSE.

!Print *, '1posit', posit, 'posit2', posit2
!Initializes counters
m = 1
t = 0.

Do
!	Print *, 'm And', m
	!Calculate distance between particles and the forces exerted on the particles
	Call force_calc(posit, fij, deltap)
!	Print *, 'acc_count', acc_count
!	Print *, '2posit', posit, 'posit2', posit2
!	posit2 = posit + vel*tstep + (fij/2/mass*(tstep**2))

	!Print *, '3posit', posit, 'posit2', posit2


	Do j1 = 1, num
		Do j2 = 1, 3
!			if(posit(j1,j2) > 10*side) then
!			Print *, 'Error - particle1'
!			Exit
!			Endif

!			if(posit2(j1,j2) > 10*side) then
!			Print *, 'Error - particle2'
!			Exit
!			Endif
		posit2(j1,j2) = posit(j1,j2) + vel(j1,j2)*tstep + (fij(j1,j2)/2/mass*(tstep**2))
!		Print *, '   posit2', posit2(j1,j2), 'posit', posit(j1,j2), 'vel', vel(j1,j2), 'fij', fij(j1,j2) !'tstep', tstep

		if(posit2(j1,j2) > side) then
			j3 = posit2(j1,j2)/side
			posit2(j1,j2) = posit2(j1,j2) - (j3 * side)
!			Print *, 'posit2', posit2(j1,j2)
!			vel(j1,j2) = - vel(j1,j2)
!			Print *, 'neg vel', vel(j1, j2)
		End if
!		elseif(posit2(j1,j2) == side) then
!			vel(j1,j2) = -vel(j1, j2)
!			Print *, 'negvel', vel(j1, j2)
		if (posit2(j1,j2) < 0.0) then
			j3 = -(posit2(j1,j2)/side) + 1
			posit2(j1,j2) = posit2(j1,j2) + (j3 * side)
!			Print *, 'neg fix', posit2(j1,j2), j3, (j3*side)
		End if

		End Do
	End Do

	!Calculates the force after the first time step, used to calculate vel2
	Call force_calc(posit2, fij2, deltap2)
!	Print *, 'acc_count', acc_count
!	vel2 = vel + ((fij+fij2)/2/mass*tstep)
!	Print *, 'vel', vel, 'vel2', vel2
	!This calculation no longer needs to be performed as the velocities
	! are in essense 'reset' below by using the subroutine init_v0

	!This calculation rescales the velocity to achieve the desired
	! temperature for the simulation, then recalculates the
	! rmsvel2 and the thermodynamic variables (temp2, press2, en2)
!	vel2 = sqrt(kbT/temp2)*vel2

	!After futher research, it appears the above method, scaling the
	! velocities, is not desired and the preferred method is to reset
	! the velocities by utilizing the same method for initializing
	! the velocities by using the Gaussian distribution as shown
	! below
	Call thermostat_v(vel2)

	Call vel_rms(vel, rmsvel)
!	Print *, 'The rms vel =', rmsvel
	Call vel_rms(vel2, rmsvel2)
!	Print *, 'The rms vel2 =', rmsvel2

	!Calculates the thermodynamic properties for the initialized variables
	Call thermo_calc(rmsvel, temp, press, en, fij, deltap, posit)
!	Print *, 'rmsvel =', rmsvel
!	Print *, 'temp =', temp
!	Print *, 'press =', press
!	Print *, 'en =', en

	!Calculates the thermodynamic properties for the initialized variables
	Call thermo_calc(rmsvel2, temp2, press2, en2, fij2, deltap2, posit2)
!	Print *, 'rmsvel2 =', rmsvel2
!	Print *, 'temp2 =', temp2
!	Print *, 'press2 =', press2
!	Print *, 'en2 =', en2


	!The following 'equilibrium' variables calculates the difference
	! between the variables: temp, press, en, and rmsvel for before
	! the time step and after the time step
!	equilT = (temp2 - temp)/temp2
!	Print *, '   The equilT =', equilT

!	equilP = (press2 - press)/press2
!	Print *, '   The equilP =', equilP

!	equilE = (en2 - en)/en2
!	Print *, '   The equilE =', equilE

!	equilVel = (rmsvel2 - rmsvel)/rmsvel2
!	Print *, '   The equilVel =', equilVel

	!The pressure is pressure is printed to compare at different densities8
!	Print *, 'Press2:', press2

	if(t >= tmax) then
		Exit
	else
		!sets vel and posit values for next time step
		vel = vel2
		posit = posit2
		!adds to the counter for the number of iterations completed
		m = m + 1
		t = t + tstep
		Cycle
	End if
End Do

End Subroutine Andersen_sim


!Subroutine rms_displace calculates the root mean square displacement
! for a given timestep
Subroutine displace_rms(rms_dis)
Use mdproj
real(Kind = DP), dimension(3):: rms_dis
integer :: loop1

rms_dis(1) = 0
rms_dis(2) = 0
rms_dis(3) = 0

Do loop1 = 1, num
!Takes the square of the displacement of each particle for the time step
rms_dis(1) = (posit2(loop1,1)-posit0(loop1,1))**2 + rms_dis(1)
rms_dis(2) = (posit2(loop1,2)-posit0(loop1,2))**2 + rms_dis(2)
rms_dis(3) = (posit2(loop1,3)-posit0(loop1,3))**2 + rms_dis(3)
End Do

!Divides the sum of the squared displacements by the number of molecules
rms_dis(1) = rms_dis(1)/num
rms_dis(2) = rms_dis(2)/num
rms_dis(3) = rms_dis(3)/num

!Print *, 'rms_dis', rms_dis
!Print *, 'posit0', posit0
!Print *, posit2

End Subroutine displace_rms


!Reset_v subroutine simulates collision with heat bath, resets a fraction
! of the velocities by the set v_frac by using the original technique
! utilized for the initialization of velocities
Subroutine thermostat_v(velin)

Use mdproj
real :: mean, stdDev, u, v, s, spare, mul, outvar, rand_num_part
integer :: dovar1, dovar2, loopcount, rand_part
real (Kind = DP), dimension(num, 3):: velin
!integer, dimension(loopcount) :: already_used
logical :: isSpareReady = .FALSE.

mean = 0.
stdDev = 1.
loopcount = num * v_frac

!Do loops simulate collision with heat bath for v_frac of particles
Do dovar1 = 1, loopcount

	Call random_number(rand_num_part)
	rand_part = num * rand_num_part

	Do dovar2 = 1, 3
	isSpareReady = .FALSE.

	!If statement handles generation of the boltmann distribution
	! for the initial veloities
	if (isSpareReady) then
		isSpareReady = .FALSE.
		outvar = spare * stdDev + mean
	else
		Do
			call random_number(u)
			u = u * 2 - 1
			call random_number(v)
			v = v * 2 - 1
			s = u*u + v*v

			if((s < 1.) .AND. (s /= 0)) then
				Exit
			End if
		End Do

		mul = sqrt(-2.0*log(s)/s)
		spare = v * mul
		isSpareREady=.TRUE.
		outvar = mean + stdDev * u * mul
	End if

	velin(rand_part,dovar2)=sqrt(kbT/mass)*outvar
!	Print *, 'part', rand_part, 'direction', dovar2, 'vel', velin(rand_part,dovar2)

	End Do
End Do

End Subroutine thermostat_v

Subroutine vel_analysis(velarrayin, velavgout, velminout, velmaxout, stddevout)
Use mdproj
real (Kind = DP), dimension(num, 3):: velarrayin
real (Kind = DP):: velavgout, velminout, velmaxout, stddevout, vel1
integer::counter_a, counter_b

velminout = 10000000
velmaxout = -10000000

Do counter_a = 1, num
	vel1 = velarrayin(counter_a, 1) + velarrayin(counter_a, 2) + velarrayin(counter_a, 3)
	if(velarrayin(counter_a, 1) < velminout) velminout = velarrayin(counter_a, 1)
	if(velarrayin(counter_a, 2) < velminout) velminout = velarrayin(counter_a, 2)
	if(velarrayin(counter_a, 3) < velminout) velminout = velarrayin(counter_a, 3)

	if(velarrayin(counter_a, 1) > velmaxout) velmaxout = velarrayin(counter_a, 1)
	if(velarrayin(counter_a, 2) > velmaxout) velmaxout = velarrayin(counter_a, 2)
	if(velarrayin(counter_a, 3) > velmaxout) velmaxout = velarrayin(counter_a, 3)

End Do

velavgout = vel1 / num / 3
vel1 = 0.

Do counter_a = 1, num
	Do counter_b = 1, 3
		vel1 = vel1 + (velarrayin(counter_a, counter_b) - velavgout)**2
	End Do
End Do

stddevout = sqrt(vel1 / num / 3)

End Subroutine vel_analysis

Subroutine Uij_calc(positin, Uij_out, deltap1)
Use mdproj
real (KIND = DP), dimension(num, 3) :: positin, Uij_out
real (Kind = DP), dimension(num, num) :: deltap1
real (Kind = DP):: u, ux, uy, uz, cx, cy, cz
integer :: count1, count2, i, j
real :: cmin, cmax
!logical :: truthtest
!infinite = HUGE(infinite)

Do count1 = 1, num
	Do count2 = 1, 3
		Uij_out(count1,count2)=0.
	End Do
End Do

cmin = 10
cmax = -10

!Print *, '0positin', positin

Do i = 1, num
	Do j = 1, num
		if(i == j) then
!			Print *, "Entered the if at i=", i, 'and j=', j
			deltap1(i,j) = 0.
!			Print *, '1positin', positin(i,1), 'j', positin(j,1)
		else
!			Print *, '1posit', posit, 'posit2', posit2
!			Print *, '2positin', positin(i,1), 'j', positin(j,1)
!			Print *, "Entered the else if at i=", i, 'and j=', j

			cx = positin(i,1) - positin(j,1)
			cy = positin(i,2) - positin(j,2)
			cz = positin(i,3) - positin(j,3)
!			Print *, 'cx', cx, 'cy', cy, 'cz', cz

!			if(cx > (side/2)) cx = cx - (side/2)
!			if(cy > (side/2)) cy = cy - (side/2)
!			if(cz > (side/2)) cz = cz - (side/2)
!			if(cx < (-side/2)) cx = cx + (side/2)
!			if(cy < (-side/2)) cy = cy + (side/2)
!			if(cz < (-side/2)) cz = cz + (side/2)

			if(cx < cmin) cmin = cx
			if(cy < cmin) cmin = cy
			if(cz < cmin) cmin = cz
			if(cx > cmax) cmax = cx
			if(cy > cmax) cmax = cy
			if(cz > cmax) cmax = cz

!			Print *, '3positin', positin(i,1), 'j', positin(j,1), 'i,j', i, j
			deltap1(i,j) = sqrt((cx)**2 + (cy)**2 + (cz)**2)
!			Print *, 'deltap1', deltap1(i,j), 'i', i, 'j', j

			if((deltap1(i,j) == 0.) .OR. (deltap1(i,j) >= cutoff)) then
!			if(deltap1(i,j) >= cutoff) then
				Uij_out(i,1) = Uij_out(i,1) + 0.0
				Uij_out(i,2) = Uij_out(i,2) + 0.0
				Uij_out(i,3) = Uij_out(i,3) + 0.0
!				Print *, 'deltap1', deltap1(i,j), 'i', i, 'j', j
!				Print *, 'cx', cx, 'cy', cy, 'cz', cz
			else
!				Print *, "the program entered the else at i=", i, 'and j=', j
				cx = cx/deltap1(i,j)
				cy = cy/deltap1(i,j)
				cz = cz/deltap1(i,j)
				u = 4/kbT*((sigma**12/(deltap1(i,j)**12)) - (sigma**6/(deltap1(i,j)**6)))
				ux = u * cx
				uy = u * cy
				uz = u * cz
				Uij_out(i,1) = Uij_out(i,1) + ux
				Uij_out(i,2) = Uij_out(i,2) + uy
				Uij_out(i,3) = Uij_out(i,3) + uz
!				Print *, 'i', i, 'j', j
!				Print *, 'ux', ux, 'uy', uy, 'uz', uz
!				Print *, 'The sum potenetial en for particle i are:', Uij_out(i,1), Uij_out(i,2), Uij_out(i,3)
!				Print *, ' '
!				Print *, ' deltap', deltap(i,j)
			End if
		End if
!		Print *, "the program completed the if statement at i=", i, 'and j=', j
	End Do
!	Print *, 'i', i, 'fij_outx', fij_out(i,1), 'fij_outy', fij_out(i,2), 'fij_outz', fij_out(i,3)
End Do
!Print *, 'cmin', cmin, 'cmax', cmax
!Print *, '2posit', posit, 'posit2', posit2

End Subroutine Uij_calc


