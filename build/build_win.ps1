$version="2.0"
$year=2026
dir .
del output\\semisim.exe
jpackage --type app-image `
	--app-version $version `
	--copyright "Brandon Li ($year)" `
	--name SemiSim `
	--input ..\\target\\ `
	--dest output\\ `
	--main-class electrodynamics.SemiSim `
	--main-jar SemiSim-$version.jar