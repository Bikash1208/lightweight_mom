$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2025
		MinorVer=2
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/11/2026 21:47:15')
			I(1, 'Host', 'AEM-LAB1')
			I(1, 'Processor', '48')
			I(1, 'OS', 'NT 10.0')
			I(1, 'Product', 'HFSS Version 2025.2.0')
		$end 'StartInfo'
		$begin 'TotalInfo'
		$end 'TotalInfo'
		GroupOptions=8
		TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Executing From\', \'C:\\\\Program Files\\\\ANSYS Inc\\\\v252\\\\AnsysEM\\\\HFSSCOMENGINE.exe\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='HPC'
			$begin 'StartInfo'
				I(1, 'Type', 'Manual')
				I(1, 'Distribution Types', 'Variations, Frequencies, Domain Solver, Transient Excitations, Direct Solver, Iterative Solver, Solver, Mesh Assembly, Mesher')
				I(1, 'MPI Vendor', 'Intel')
				I(1, 'MPI Version', '2021')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(0, ' ')
			$end 'TotalInfo'
			GroupOptions=0
			TaskDataOptions(Memory=8)
			ProfileItem('Two Level', 0, 0, 0, 0, 0, 'I(1, 0, \'Disabled\')', false, true)
			ProfileItem('Machine', 0, 0, 0, 0, 0, 'I(6, 1, \'Name\', \'AEM-LAB1\', 1, \'Memory\', \'256 GB\', 3, \'RAM Limit\', 90, \'%f%%\', 2, \'Tasks\', 1, false, 2, \'Cores\', 5, false, 1, \'Free Disk Space\', \'311 GB\')', false, true)
		$end 'ProfileGroup'
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Allow off core\', \'True\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Solution Basis Order\', \'1\')', false, true)
		ProfileItem('Design Validation', 0, 0, 0, 0, 0, 'I(1, 0, \'Elapsed time : 00:00:00 , HFSS ComEngine Memory : 189 M\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'Perform full validations with standard port validations\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Initial Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 21:47:18')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:18:10')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 601, 0, 810, 0, 5460000, 'I(3, 2, \'Tetrahedra\', 4254476, false, 1, \'Type\', \'TAU\', 2, \'Cores\', 5, false)', true, true)
			ProfileItem('Coarsen', 355, 0, 355, 0, 5460000, 'I(1, 2, \'Tetrahedra\', 1178785, false)', true, true)
			ProfileItem('Lambda Refine', 26, 0, 29, 0, 1147940, 'I(1, 2, \'Tetrahedra\', 1200262, false)', true, true)
			ProfileItem('Simulation Setup', 23, 0, 23, 0, 2811864, 'I(2, 2, \'Tetrahedra\', 788278, false, 1, \'Disk\', \'0 Bytes\')', true, true)
			ProfileItem('Port Adapt', 14, 0, 14, 0, 2830132, 'I(2, 2, \'Tetrahedra\', 788278, false, 1, \'Disk\', \'76.1 KB\')', true, true)
			ProfileItem('Port Refine', 58, 0, 61, 0, 1500200, 'I(1, 2, \'Tetrahedra\', 1196263, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 22:05:28')
			$end 'StartInfo'
			$begin 'TotalInfo'
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			$begin 'ProfileGroup'
				MajorVer=2025
				MinorVer=2
				Name='Adaptive Pass 1'
				$begin 'StartInfo'
					I(1, 'Frequency', '5.8GHz')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 21, 0, 22, 0, 2805656, 'I(2, 2, \'Tetrahedra\', 789019, false, 1, \'Disk\', \'185 KB\')', true, true)
				ProfileItem('Matrix Assembly', 90, 0, 162, 0, 7711216, 'I(3, 2, \'Tetrahedra\', 789019, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'305 MB\')', true, true)
				ProfileItem('Matrix Solve', 164, 0, 549, 0, 16871940, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 4114809, false, 3, \'Matrix bandwidth\', 17.455, \'%5.1f\', 1, \'Disk\', \'15.7 MB\')', true, true)
				ProfileItem('Field Recovery', 18, 0, 37, 0, 16871940, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'99.4 MB\')', true, true)
				ProfileItem('Data Transfer', 4, 0, 4, 0, 208860, 'I(1, 0, \'Adaptive Pass 1\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415207')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'4114809\', 1, \'Matrixbw\', \'17.456499\', 1, \'Matrixnnz\', \'71830145\', 1, \'Rootdim\', \'5740\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'365883\', 1, \'Factornnz\', \'1103671343\', 1, \'Factorestflops\', \'4516490000000\', 1, \'Fbsestflops\', \'3651260000\', 1, \'Rootfactestflops\', \'63039800000\', 1, \'Rootfbsestflops\', \'16473800\', 1, \'Analysistimesec\', \'60.868400\', 1, \'Analysismemkb\', \'2106840.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'91.395200\', 1, \'Factorizationmentotalkb\', \'13414600.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'8.605700\', 1, \'Fbscputimesec\', \'8.605700\', 1, \'Fbsmemorytotalkb\', \'11029600.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'13414620\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
			$end 'ProfileGroup'
			ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
			$begin 'ProfileGroup'
				MajorVer=2025
				MinorVer=2
				Name='Adaptive Pass 2'
				$begin 'StartInfo'
					I(1, 'Frequency', '5.8GHz')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 94, 0, 99, 0, 1799824, 'I(1, 2, \'Tetrahedra\', 1432974, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 23, 0, 24, 0, 3370608, 'I(2, 2, \'Tetrahedra\', 1014800, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 189, 0, 311, 0, 10283960, 'I(3, 2, \'Tetrahedra\', 1014800, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'430 MB\')', true, true)
				ProfileItem('Matrix Solve', 331, 0, 1227, 0, 26837032, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 5519949, false, 3, \'Matrix bandwidth\', 18.4411, \'%5.1f\', 1, \'Disk\', \'21.1 MB\')', true, true)
				ProfileItem('Field Recovery', 34, 0, 58, 0, 26837032, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'45 MB\')', true, true)
				ProfileItem('Data Transfer', 31, 0, 8, 0, 208976, 'I(1, 0, \'Adaptive Pass 2\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415875')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'5519949\', 1, \'Matrixbw\', \'18.443800\', 1, \'Matrixnnz\', \'101808899\', 1, \'Rootdim\', \'8250\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'486713\', 1, \'Factornnz\', \'1909663529\', 1, \'Factorestflops\', \'11569100000000\', 1, \'Fbsestflops\', \'6410740000\', 1, \'Rootfactestflops\', \'187172000000\', 1, \'Rootfbsestflops\', \'34031250\', 1, \'Analysistimesec\', \'89.623800\', 1, \'Analysismemkb\', \'2975260.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'216.014000\', 1, \'Factorizationmentotalkb\', \'24126200.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'18.434500\', 1, \'Fbscputimesec\', \'18.434500\', 1, \'Fbsmemorytotalkb\', \'19209100.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'24126190\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.14694, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
