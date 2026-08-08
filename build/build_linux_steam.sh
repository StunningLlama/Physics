version=$(cat version.txt)
year=$(cat year.txt)
rm -rf output_linux
jpackage --type app-image \
	--app-version $version \
	--copyright "Brandon Li ($year)" \
	--description "Brandon's Semiconductor Simulator" \
	--name SemiSim \
	--icon ../images/icon.png \
	--vendor "Brandon Li" \
	--input ../target/ \
	--dest output_linux \
	--main-class electrodynamics.SemiSim \
	--main-jar SemiSim-$version.jar \
	--java-options -XX:-TieredCompilation \
	--java-options -XX:CompileThresholdScaling=0.25

sh steam/ContentBuilder/scripts/build_linux.sh