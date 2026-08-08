$version = Get-Content ./version.txt -Raw
$year = Get-Content ./year.txt -Raw
del output_win_installer -Recurse -Force
jpackage --type exe `
	--app-version $version `
	--copyright "Brandon Li ($year)" `
	--description "Brandon's Semiconductor Simulator" `
	--name SemiSim `
	--icon ..\\images\\icon.ico `
	--vendor "Brandon Li" `
	--input ..\\target\\ `
	--dest output_win_installer\\ `
	--main-class electrodynamics.SemiSim `
	--main-jar SemiSim-$version.jar `
	--java-options -XX:-TieredCompilation `
	--java-options -XX:CompileThresholdScaling=0.25 `
	--about-url "https://brandonli.net/semisim" `
	--file-associations semisim.properties `
	--license-file ..\\license.txt `
	--win-dir-chooser `
	--win-menu `
	--win-menu-group "SemiSim" `
	--win-shortcut `
	--win-shortcut-prompt