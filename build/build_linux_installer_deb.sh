version=$(cat version.txt)
year=$(cat year.txt)
rm -rf output_linux_installer_deb
jpackage --type deb \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--description "Brandon's Semiconductor Simulator" \
	--name SemiSim \
	--icon ../images/icon.png \
	--vendor "Brandon Li" \
	--input ../target/ \
	--dest output_linux_installer_deb \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar \
	--java-options -XX:-TieredCompilation \
	--java-options -XX:CompileThresholdScaling=0.25 \
	--about-url "https://brandonli.net/semisim" \
	--file-associations semisim.properties \
	--license-file ../license.txt \
	--linux-shortcut