$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2025
		MinorVer=2
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/11/2026 21:47:14')
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
				I(1, 'Time', '11/11/2026 21:47:16')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:17:58')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 594, 0, 784, 0, 5464000, 'I(3, 2, \'Tetrahedra\', 4254476, false, 1, \'Type\', \'TAU\', 2, \'Cores\', 5, false)', true, true)
			ProfileItem('Coarsen', 356, 0, 356, 0, 5464000, 'I(1, 2, \'Tetrahedra\', 1178785, false)', true, true)
			ProfileItem('Lambda Refine', 30, 0, 34, 0, 1153144, 'I(1, 2, \'Tetrahedra\', 1200552, false)', true, true)
			ProfileItem('Simulation Setup', 21, 0, 21, 0, 2812672, 'I(2, 2, \'Tetrahedra\', 788576, false, 1, \'Disk\', \'0 Bytes\')', true, true)
			ProfileItem('Port Adapt', 13, 0, 13, 0, 2831004, 'I(2, 2, \'Tetrahedra\', 788576, false, 1, \'Disk\', \'76.1 KB\')', true, true)
			ProfileItem('Port Refine', 52, 0, 55, 0, 1491308, 'I(1, 2, \'Tetrahedra\', 1196578, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 22:05:15')
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
				ProfileItem('Simulation Setup ', 22, 0, 22, 0, 2805100, 'I(2, 2, \'Tetrahedra\', 789334, false, 1, \'Disk\', \'185 KB\')', true, true)
				ProfileItem('Matrix Assembly', 91, 0, 161, 0, 7736824, 'I(3, 2, \'Tetrahedra\', 789334, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'306 MB\')', true, true)
				ProfileItem('Matrix Solve', 164, 0, 555, 0, 17021128, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 4116781, false, 3, \'Matrix bandwidth\', 17.4569, \'%5.1f\', 1, \'Disk\', \'15.7 MB\')', true, true)
				ProfileItem('Field Recovery', 19, 0, 38, 0, 17021128, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'99.4 MB\')', true, true)
				ProfileItem('Data Transfer', 3, 0, 3, 0, 207660, 'I(1, 0, \'Adaptive Pass 1\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415195')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'4116781\', 1, \'Matrixbw\', \'17.458401\', 1, \'Matrixnnz\', \'71872549\', 1, \'Rootdim\', \'5740\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'366073\', 1, \'Factornnz\', \'1114851227\', 1, \'Factorestflops\', \'4635530000000\', 1, \'Fbsestflops\', \'3697870000\', 1, \'Rootfactestflops\', \'63039800000\', 1, \'Rootfbsestflops\', \'16473800\', 1, \'Analysistimesec\', \'59.200400\', 1, \'Analysismemkb\', \'2108040.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'91.243200\', 1, \'Factorizationmentotalkb\', \'13819300.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'9.172560\', 1, \'Fbscputimesec\', \'9.172560\', 1, \'Fbsmemorytotalkb\', \'11154700.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'13819256\', 1, \'Maxdiskkb\', \'0\')', false, true)
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
				ProfileItem('Adaptive Refine', 108, 0, 114, 0, 1803600, 'I(1, 2, \'Tetrahedra\', 1433386, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 34, 0, 35, 0, 3371316, 'I(2, 2, \'Tetrahedra\', 1014932, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 191, 0, 314, 0, 10285916, 'I(3, 2, \'Tetrahedra\', 1014932, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'430 MB\')', true, true)
				ProfileItem('Matrix Solve', 324, 0, 1189, 0, 26286272, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 5, false, 2, \'Matrix size\', 5520633, false, 3, \'Matrix bandwidth\', 18.4346, \'%5.1f\', 1, \'Disk\', \'21.1 MB\')', true, true)
				ProfileItem('Field Recovery', 34, 0, 52, 0, 26286272, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'45 MB\')', true, true)
				ProfileItem('Data Transfer', 32, 0, 10, 0, 207820, 'I(1, 0, \'Adaptive Pass 2\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415882')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'5520633\', 1, \'Matrixbw\', \'18.437401\', 1, \'Matrixnnz\', \'101785929\', 1, \'Rootdim\', \'8358\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'487568\', 1, \'Factornnz\', \'1875598867\', 1, \'Factorestflops\', \'11067200000000\', 1, \'Fbsestflops\', \'6303840000\', 1, \'Rootfactestflops\', \'194619000000\', 1, \'Rootfbsestflops\', \'34928082\', 1, \'Analysistimesec\', \'95.253600\', 1, \'Analysismemkb\', \'2974660.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'5\', 1, \'Factorizationtimesec\', \'207.747000\', 1, \'Factorizationmentotalkb\', \'23601000.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'5\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'14.287900\', 1, \'Fbscputimesec\', \'14.287900\', 1, \'Fbsmemorytotalkb\', \'18656000.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'23600957\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.175275, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
