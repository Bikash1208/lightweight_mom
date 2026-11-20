$begin 'Profile'
	$begin 'ProfileGroup'
		MajorVer=2026
		MinorVer=1
		Name='Solution Process'
		$begin 'StartInfo'
			I(1, 'Start Time', '11/20/2026 18:24:26')
			I(1, 'Host', 'AEM-LAB1')
			I(1, 'Processor', '48')
			I(1, 'OS', 'NT 10.0')
			I(1, 'Product', 'HFSS Version 2026.1.0')
		$end 'StartInfo'
		$begin 'TotalInfo'
		$end 'TotalInfo'
		GroupOptions=8
		TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Executing From\', \'C:\\\\Program Files\\\\ANSYS Inc\\\\v261\\\\AnsysEM\\\\HFSSCOMENGINE.exe\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2026
			MinorVer=1
			Name='HPC'
			$begin 'StartInfo'
				I(1, 'Type', 'Auto')
				I(1, 'MPI Vendor', 'Intel')
				I(1, 'MPI Version', '2021')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(0, ' ')
			$end 'TotalInfo'
			GroupOptions=0
			TaskDataOptions(Memory=8)
			ProfileItem('Machine', 0, 0, 0, 0, 0, 'I(5, 1, \'Name\', \'AEM-LAB1\', 1, \'Memory\', \'256 GB\', 3, \'RAM Limit\', 89.9998, \'%f%%\', 2, \'Cores\', 44, false, 1, \'Free Disk Space\', \'304 GB\')', false, true)
		$end 'ProfileGroup'
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Allow off core\', \'True\')', false, true)
		ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 1, \'Solution Basis Order\', \'1\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2026
			MinorVer=1
			Name='Design Validation'
			$begin 'StartInfo'
				I(1, 'Level', 'Perform full validations with standard port validations')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:00:00')
				I(1, 'Memory', '107 M')
			$end 'TotalInfo'
			GroupOptions=0
			TaskDataOptions(Memory=8)
		$end 'ProfileGroup'
		ProfileItem('Mesh Fusion Preprocessing', 0, 0, 0, 0, 0, 'I(1, 0, \'Elapsed time : 00:00:00\')', false, true)
		$begin 'ProfileGroup'
			MajorVer=2026
			MinorVer=1
			Name='Initial Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/20/2026 18:24:27')
			$end 'StartInfo'
			$begin 'TotalInfo'
				I(1, 'Elapsed Time', '00:00:03')
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			ProfileItem('Mesh', 0, 0, 0, 0, 44916, 'I(2, 1, \'Type\', \'Phi\', 2, \'Tetrahedra\', 2824, false)', true, true)
			ProfileItem('Post', 1, 0, 1, 0, 50660, 'I(1, 2, \'Tetrahedra\', 3983, false)', true, true)
			ProfileItem('Lambda Refine', 1, 0, 1, 0, 34780, 'I(1, 2, \'Tetrahedra\', 9572, false)', true, true)
		$end 'ProfileGroup'
		$begin 'ProfileGroup'
			MajorVer=2026
			MinorVer=1
			Name='Adaptive Meshing'
			$begin 'StartInfo'
				I(1, 'Time', '11/20/2026 18:24:30')
			$end 'StartInfo'
			$begin 'TotalInfo'
			$end 'TotalInfo'
			GroupOptions=4
			TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 1'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 210700, 'I(2, 2, \'Tetrahedra\', 8336, false, 1, \'Disk\', \'3.36 KB\')', true, true)
				ProfileItem('Matrix Assembly', 1, 0, 10, 0, 370756, 'I(2, 2, \'Tetrahedra\', 8336, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 154, 0, 1809, 0, 1.00708e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 49228, false, 3, \'Matrix bandwidth\', 41.0177, \'%5.1f\', 2, \'Eigen iterations\', 179, false, 1, \'Disk\', \'15 MB\')', true, true)
				ProfileItem('Field Recovery', 1, 0, 49, 0, 1.00708e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'80 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
			$end 'ProfileGroup'
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 2'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 0, 0, 0, 0, 36776, 'I(1, 2, \'Tetrahedra\', 12056, false)', true, true)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 218344, 'I(2, 2, \'Tetrahedra\', 10688, false, 1, \'Disk\', \'2.98 KB\')', true, true)
				ProfileItem('Matrix Assembly', 1, 0, 10, 0, 424232, 'I(2, 2, \'Tetrahedra\', 10688, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 103, 0, 1184, 0, 1.28096e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 63828, false, 3, \'Matrix bandwidth\', 41.5064, \'%5.1f\', 2, \'Eigen iterations\', 109, false, 1, \'Disk\', \'4.46 MB\')', true, true)
				ProfileItem('Field Recovery', 1, 0, 60, 0, 1.28096e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'22.7 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileFootnote('I(1, 3, \'Max Delta Freq. %\', 4.38747, \'%.5f\')', 0)
			$end 'ProfileGroup'
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 3'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 0, 0, 0, 0, 40924, 'I(1, 2, \'Tetrahedra\', 15267, false)', true, true)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 226364, 'I(2, 2, \'Tetrahedra\', 13808, false, 1, \'Disk\', \'7.19 KB\')', true, true)
				ProfileItem('Matrix Assembly', 2, 0, 10, 0, 491336, 'I(2, 2, \'Tetrahedra\', 13808, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 150, 0, 1653, 0, 1.64237e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 83122, false, 3, \'Matrix bandwidth\', 41.8875, \'%5.1f\', 2, \'Eigen iterations\', 116, false, 1, \'Disk\', \'5.89 MB\')', true, true)
				ProfileItem('Field Recovery', 1, 0, 79, 0, 1.64237e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'30.1 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileFootnote('I(1, 3, \'Max Delta Freq. %\', 5.04552, \'%.5f\')', 0)
			$end 'ProfileGroup'
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 4'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 0, 0, 1, 0, 45432, 'I(1, 2, \'Tetrahedra\', 19427, false)', true, true)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 237816, 'I(2, 2, \'Tetrahedra\', 17828, false, 1, \'Disk\', \'7.19 KB\')', true, true)
				ProfileItem('Matrix Assembly', 2, 0, 11, 0, 579836, 'I(2, 2, \'Tetrahedra\', 17828, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 646, 0, 6757, 0, 2.16554e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 107964, false, 3, \'Matrix bandwidth\', 42.1759, \'%5.1f\', 2, \'Eigen iterations\', 83, false, 1, \'Disk\', \'7.58 MB\')', true, true)
				ProfileItem('Field Recovery', 2, 0, 99, 0, 2.16554e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'38.8 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileFootnote('I(1, 3, \'Max Delta Freq. %\', 2.18998, \'%.5f\')', 0)
			$end 'ProfileGroup'
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 5'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 0, 0, 0, 0, 47788, 'I(1, 2, \'Tetrahedra\', 21956, false)', true, true)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 242860, 'I(2, 2, \'Tetrahedra\', 20251, false, 1, \'Disk\', \'9.49 KB\')', true, true)
				ProfileItem('Matrix Assembly', 3, 0, 11, 0, 631288, 'I(2, 2, \'Tetrahedra\', 20251, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 342, 0, 3220, 0, 2.40628e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 122790, false, 3, \'Matrix bandwidth\', 42.2278, \'%5.1f\', 2, \'Eigen iterations\', 157, false, 1, \'Disk\', \'4.52 MB\')', true, true)
				ProfileItem('Field Recovery', 2, 0, 110, 0, 2.40628e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'23.6 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileFootnote('I(1, 3, \'Max Delta Freq. %\', 1.20786, \'%.5f\')', 0)
			$end 'ProfileGroup'
			$begin 'ProfileGroup'
				MajorVer=2026
				MinorVer=1
				Name='Adaptive Pass 6'
				$begin 'StartInfo'
					I(0, 'Eigenmode Solution')
				$end 'StartInfo'
				$begin 'TotalInfo'
					I(0, ' ')
				$end 'TotalInfo'
				GroupOptions=0
				TaskDataOptions('CPU Time'=8, Memory=8, 'Real Time'=8)
				ProfileItem('Adaptive Refine', 1, 0, 1, 0, 53396, 'I(1, 2, \'Tetrahedra\', 27080, false)', true, true)
				ProfileItem('Simulation Setup', 0, 0, 0, 0, 255880, 'I(2, 2, \'Tetrahedra\', 25114, false, 1, \'Disk\', \'9.49 KB\')', true, true)
				ProfileItem('Matrix Assembly', 4, 0, 13, 0, 735492, 'I(2, 2, \'Tetrahedra\', 25114, false, 1, \'Disk\', \'0 Bytes\')', true, true)
				ProfileItem('Matrix Solve', 963, 0, 9414, 0, 3.05476e+06, 'I(6, 1, \'Type\', \'DCA\', 2, \'Cores\', 44, false, 2, \'Matrix size\', 152728, false, 3, \'Matrix bandwidth\', 42.3492, \'%5.1f\', 2, \'Eigen iterations\', 83, false, 1, \'Disk\', \'9.14 MB\')', true, true)
				ProfileItem('Field Recovery', 3, 0, 134, 0, 3.05476e+06, 'I(2, 2, \'Computed eigenmodes\', 20, false, 1, \'Disk\', \'47 MB\')', true, true)
				ProfileItem('', 0, 0, 0, 0, 0, 'I(1, 0, \'\')', false, true)
				ProfileFootnote('I(1, 3, \'Max Delta Freq. %\', 1.10471, \'%.5f\')', 0)
			$end 'ProfileGroup'
		$end 'ProfileGroup'
	$end 'ProfileGroup'
$end 'Profile'
