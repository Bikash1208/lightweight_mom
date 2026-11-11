$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2025
		MinorVer=2
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/11/2026 21:45:36')
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
			ProfileItem('Machine', 0, 0, 0, 0, 0, 'I(6, 1, \'Name\', \'AEM-LAB1\', 1, \'Memory\', \'256 GB\', 3, \'RAM Limit\', 90, \'%f%%\', 2, \'Tasks\', 1, false, 2, \'Cores\', 12, false, 1, \'Free Disk Space\', \'311 GB\')', false, true)
		$end 'ProfileGroup'
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Allow off core\', \'True\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Solution Basis Order\', \'1\')', false, true)
		ProfileItem('Design Validation', 0, 0, 0, 0, 0, 'I(1, 0, \'Elapsed time : 00:00:00 , HFSS ComEngine Memory : 231 M\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'Perform full validations with standard port validations\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Initial Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 21:45:38')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:17:41')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 577, 0, 967, 0, 5264000, 'I(3, 2, \'Tetrahedra\', 4183852, false, 1, \'Type\', \'TAU\', 2, \'Cores\', 12, false)', true, true)
			ProfileItem('Coarsen', 361, 0, 361, 0, 5264000, 'I(1, 2, \'Tetrahedra\', 1153984, false)', true, true)
			ProfileItem('Lambda Refine', 23, 0, 27, 0, 1127816, 'I(1, 2, \'Tetrahedra\', 1175399, false)', true, true)
			ProfileItem('Simulation Setup', 23, 0, 23, 0, 2765088, 'I(2, 2, \'Tetrahedra\', 774586, false, 1, \'Disk\', \'0 Bytes\')', true, true)
			ProfileItem('Port Adapt', 15, 0, 15, 0, 2779592, 'I(2, 2, \'Tetrahedra\', 774586, false, 1, \'Disk\', \'76.6 KB\')', true, true)
			ProfileItem('Port Refine', 49, 0, 53, 0, 1491900, 'I(1, 2, \'Tetrahedra\', 1170663, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2025
			MinorVer=2
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/11/2026 22:03:20')
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
				ProfileItem('Simulation Setup ', 19, 0, 21, 0, 2756336, 'I(2, 2, \'Tetrahedra\', 774732, false, 1, \'Disk\', \'185 KB\')', true, true)
				ProfileItem('Matrix Assembly', 84, 0, 155, 0, 7585216, 'I(3, 2, \'Tetrahedra\', 774732, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'301 MB\')', true, true)
				ProfileItem('Matrix Solve', 101, 0, 586, 0, 17344844, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 12, false, 2, \'Matrix size\', 4044339, false, 3, \'Matrix bandwidth\', 17.4797, \'%5.1f\', 1, \'Disk\', \'15.4 MB\')', true, true)
				ProfileItem('Field Recovery', 17, 0, 45, 0, 17344844, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'97.6 MB\')', true, true)
				ProfileItem('Data Transfer', 2, 0, 2, 0, 251328, 'I(1, 0, \'Adaptive Pass 1\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415007')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'4044339\', 1, \'Matrixbw\', \'17.481300\', 1, \'Matrixnnz\', \'70700202\', 1, \'Rootdim\', \'5512\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'359280\', 1, \'Factornnz\', \'1101023478\', 1, \'Factorestflops\', \'4612430000000\', 1, \'Fbsestflops\', \'3641740000\', 1, \'Rootfactestflops\', \'55822200000\', 1, \'Rootfbsestflops\', \'15191072\', 1, \'Analysistimesec\', \'43.605200\', 1, \'Analysismemkb\', \'2073510.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'12\', 1, \'Factorizationtimesec\', \'46.192300\', 1, \'Factorizationmentotalkb\', \'15846000.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'12\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'8.187010\', 1, \'Fbscputimesec\', \'8.187010\', 1, \'Fbsmemorytotalkb\', \'11599300.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'15845969\', 1, \'Maxdiskkb\', \'0\')', false, true)
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
				ProfileItem('Adaptive Refine', 92, 0, 102, 0, 1785488, 'I(1, 2, \'Tetrahedra\', 1403086, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 45, 0, 47, 0, 3314032, 'I(2, 2, \'Tetrahedra\', 996194, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 142, 0, 255, 0, 10108892, 'I(3, 2, \'Tetrahedra\', 996194, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'423 MB\')', true, true)
				ProfileItem('Matrix Solve', 236, 0, 1557, 0, 27291928, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 12, false, 2, \'Matrix size\', 5422539, false, 3, \'Matrix bandwidth\', 18.4507, \'%5.1f\', 1, \'Disk\', \'20.7 MB\')', true, true)
				ProfileItem('Field Recovery', 27, 0, 60, 0, 27291928, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'44.1 MB\')', true, true)
				ProfileItem('Data Transfer', 25, 0, 6, 0, 251516, 'I(1, 0, \'Adaptive Pass 2\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794415552')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'5422539\', 1, \'Matrixbw\', \'18.453400\', 1, \'Matrixnnz\', \'100064230\', 1, \'Rootdim\', \'8240\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'478649\', 1, \'Factornnz\', \'1863129940\', 1, \'Factorestflops\', \'11001400000000\', 1, \'Fbsestflops\', \'6254000000\', 1, \'Rootfactestflops\', \'186492000000\', 1, \'Rootfbsestflops\', \'33948800\', 1, \'Analysistimesec\', \'89.409700\', 1, \'Analysismemkb\', \'2924360.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'12\', 1, \'Factorizationtimesec\', \'125.547000\', 1, \'Factorizationmentotalkb\', \'28536300.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'12\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'14.129200\', 1, \'Fbscputimesec\', \'14.129200\', 1, \'Fbsmemorytotalkb\', \'19793700.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'28536325\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.229975, \'%.5f\')', 0)
			$end 'ProfileGroup'
			ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
			$begin 'ProfileGroup'
				MajorVer=2025
				MinorVer=2
				Name='Adaptive Pass 3'
				$begin 'StartInfo'
					I(1, 'Frequency', '5.8GHz')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 129, 0, 140, 0, 2046980, 'I(1, 2, \'Tetrahedra\', 1701949, false)', true, true)
				ProfileItem(' ', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileItem('Simulation Setup ', 55, 0, 56, 0, 3979224, 'I(2, 2, \'Tetrahedra\', 1280176, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Assembly', 245, 0, 387, 0, 13312940, 'I(3, 2, \'Tetrahedra\', 1280176, false, 2, \'Lumped ports\', 1, false, 1, \'Disk (output)\', \'582 MB\')', true, true)
				ProfileItem('Matrix Solve', 329, 0, 2264, 0, 39632876, 'I(5, 1, \'Type\', \'MCS\', 2, \'Cores\', 12, false, 2, \'Matrix size\', 7186929, false, 3, \'Matrix bandwidth\', 19.2108, \'%5.1f\', 1, \'Disk\', \'27.4 MB\')', true, true)
				ProfileItem('Field Recovery', 32, 0, 68, 0, 39632876, 'I(2, 2, \'Excitations\', 1, false, 1, \'Disk\', \'56.8 MB\')', true, true)
				ProfileItem('Data Transfer', 30, 0, 8, 0, 252036, 'I(1, 0, \'Adaptive Pass 3\')', true, true)
				$begin 'ProfileGroup'
					MajorVer=2025
					MinorVer=2
					Name='APIPms'
					$begin 'StartInfo'
						I(1, 'Timesinceepock', '1794416377')
					$end 'StartInfo'
					$begin 'TotalInfo'
						I(0, ' ')
					$end 'TotalInfo'
					GroupOptions=16
					TaskDataOptions(Memory=8)
					ProfileItem('solverinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Solvertype\', \'shared_memory\', 1, \'Precision\', \'mixed\', 1, \'Solversymmetry\', \'complex_sym\', 1, \'Matrixdim\', \'7186929\', 1, \'Matrixbw\', \'19.214100\', 1, \'Matrixnnz\', \'138090215\', 1, \'Rootdim\', \'10052\', 1, \'Mathtype\', \'mkl\', 1, \'Mpitasks\', \'1\', 1, \'Threadspertasks\', \'0\')', false, true)
					ProfileItem('sysinfo', 0, 0, 0, 0, 0, 'I(12, 1, \'Os\', \'win\', 1, \'Cpuid\', \'Intel(R) Xeon(R) Silver 4310 CPU @ 2.10GHz\', 1, \'CpuPhysicCores\', \'24\', 1, \'CpuLogicCores\', \'48\', 1, \'Cpufreqkhz\', \'2100000000.000000\', 1, \'Cpucachelinesizebytes\', \'64\', 1, \'Cpuestlastlevelcachesizemb\', \'30.000000\', 1, \'Cpuestgflops\', \'1344.000000\', 1, \'Memorybwestkbps\', \'170.000000\', 1, \'Numanodes\', \'2\', 1, \'Virtualmemkb\', \'137439000000.000000\', 1, \'Pagesizekb\', \'4096\')', false, true)
					ProfileItem('analysisinfo', 0, 0, 0, 0, 0, 'I(9, 1, \'Analysisstatus\', \'valid\', 1, \'Numsupernodes\', \'630130\', 1, \'Factornnz\', \'2860290000\', 1, \'Factorestflops\', \'21319800000000\', 1, \'Fbsestflops\', \'9693800000\', 1, \'Rootfactestflops\', \'338560000000\', 1, \'Rootfbsestflops\', \'50521352\', 1, \'Analysistimesec\', \'127.788000\', 1, \'Analysismemkb\', \'4025420.000000\')', false, true)
					ProfileItem('factorinfo', 0, 0, 0, 0, 0, 'I(4, 1, \'Fatorizationstatus\', \'valid\', 1, \'Factorizationnumcores\', \'12\', 1, \'Factorizationtimesec\', \'177.495000\', 1, \'Factorizationmentotalkb\', \'42500800.000000\')', false, true)
					ProfileItem('fbsinfo', 0, 0, 0, 0, 0, 'I(10, 1, \'Fbstatus\', \'valid\', 1, \'Fbstype\', \'fullsolve\', 1, \'Fbsmt\', \'false\', 1, \'Fbsmrhs\', \'false\', 1, \'Fbsnumcores\', \'12\', 1, \'Fbsnumsolvestotal\', \'1\', 1, \'Fbsnumsolves\', \'1\', 1, \'Fbsavgsolvetime1solvesec\', \'15.818400\', 1, \'Fbscputimesec\', \'15.818400\', 1, \'Fbsmemorytotalkb\', \'29929600.000000\')', false, true)
					ProfileItem('solverprofile', 0, 0, 0, 0, 0, 'I(2, 1, \'Maxmemkb\', \'42500762\', 1, \'Maxdiskkb\', \'0\')', false, true)
				$end 'ProfileGroup'
				ProfileFootnote('I(1, 3, \'Max Mag. Delta S\', 0.152284, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
