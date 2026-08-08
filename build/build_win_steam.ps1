$version = Get-Content ./version.txt -Raw
$year = Get-Content ./year.txt -Raw
del output_win -Recurse -Force
jpackage --type app-image `
	--app-version $version `
	--copyright "Brandon Li ($year)" `
	--description "Brandon's Semiconductor Simulator" `
	--name SemiSim `
	--icon ../images/icon.ico `
	--vendor "Brandon Li" `
	--input ..\\target\\ `
	--dest output_win\\ `
	--main-class electrodynamics.SemiSim `
	--main-jar SemiSim-$version.jar `
	--java-options -XX:-TieredCompilation `
	--java-options -XX:CompileThresholdScaling=0.25

steam\\ContentBuilder\\scripts\\build_win.ps1