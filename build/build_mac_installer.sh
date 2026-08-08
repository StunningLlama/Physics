version=$(cat version.txt)
year=$(cat year.txt)
rm -rf output_mac_installer/
jpackage --type dmg \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--description "Brandon's Semiconductor Simulator" \
	--name SemiSim \
	--icon ../images/icon.icns \
	--vendor "Brandon Li" \
	--input ../target/ \
	--dest output_mac_installer \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar \
	--java-options -XX:-TieredCompilation \
	--java-options -XX:CompileThresholdScaling=0.25 \
	--about-url "https://brandonli.net/semisim" \
	--file-associations semisim.properties \
	--license-file ../license.txt