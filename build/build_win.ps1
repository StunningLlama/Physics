$version="2.0"
$year=2026
dir .
del output_win -Recurse -Force
jpackage --type app-image `
	--app-version $version `
	--copyright "Brandon Li ($year)" `
	--name SemiSim `
	--input ..\\target\\ `
	--dest output_win\\ `
	--main-class electrodynamics.SemiSim `
	--main-jar SemiSim-$version.jar

steam\\ContentBuilder\\scripts\\build_win.ps1