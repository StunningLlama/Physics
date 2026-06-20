// Source: https://stackoverflow.com/questions/38931111/how-to-make-plots-in-java-like-in-matlab-same-syntax

package electrodynamics.plot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Stroke;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import electrodynamics.SemiSim;


public class MatlabChart {

    public JFreeChart chart;
    public LegendTitle legend;
    public ArrayList<Color> colors;
    public ArrayList<Stroke> strokes;
    public XYSeriesCollection dataset;

    public MatlabChart() {
        colors = new ArrayList<>();
        strokes = new ArrayList<>();
        dataset = new XYSeriesCollection();
    }

    public void plot(double[] x, double[] y, String spec, float lineWidth, String title) {
        final XYSeries series = new XYSeries(title);
        for (int i = 0; i < x.length; i++)
            series.add(x[i],y[i]);
        dataset.addSeries(series);
        FindColor(spec,lineWidth);
    }

    public XYSeries plot(String spec, float lineWidth, String title) {
        XYSeries series = new XYSeries(title);
        dataset.addSeries(series);
        FindColor(spec,lineWidth);
        return series;
    }

    public void RenderPlot() {
        // Create chart
        JFreeChart chart = null;
        if (dataset != null && dataset.getSeriesCount() > 0)
            chart = ChartFactory.createXYLineChart(null,null,null,dataset,PlotOrientation.VERTICAL,true, false, false);
        else
            System.out.println(" [!] First create a chart and add data to it. The plot is empty now!");
        // Add customization options to chart
        XYPlot plot = chart.getXYPlot();
        for (int i = 0; i < colors.size(); i++) {
            plot.getRenderer().setSeriesPaint(i, colors.get(i));
            plot.getRenderer().setSeriesStroke(i, strokes.get(i));
        }
        ((NumberAxis)plot.getDomainAxis()).setAutoRangeIncludesZero(false);
        ((NumberAxis)plot.getRangeAxis()).setAutoRangeIncludesZero(false);
        plot.setBackgroundPaint(Color.WHITE);
        legend = chart.getLegend();
        chart.removeLegend();
        this.chart = chart;
    }

    public void CheckExists() {
        if (chart == null) {
            throw new IllegalArgumentException("First plot something in the chart before you modify it.");
        }
    }

    public void font(Font font) {
        CheckExists();
        chart.getTitle().setFont(font);
        chart.getXYPlot().getDomainAxis().setLabelFont(font);
        chart.getXYPlot().getDomainAxis().setTickLabelFont(font);
        chart.getXYPlot().getRangeAxis().setLabelFont(font);
        chart.getXYPlot().getRangeAxis().setTickLabelFont(font);
        legend.setItemFont(font);
    }

    public void title(String title) {
        CheckExists();
        chart.setTitle(title);
    }

    public void xlim(double l, double u) {
        CheckExists();
        chart.getXYPlot().getDomainAxis().setRange(l, u);
    }

    public void ylim(double l, double u) {
        CheckExists();
        chart.getXYPlot().getRangeAxis().setRange(l, u);
    }

    public void xlabel(String label) {
        CheckExists();
        chart.getXYPlot().getDomainAxis().setLabel(label);
    }

    public void ylabel(String label) {
        CheckExists();
        chart.getXYPlot().getRangeAxis().setLabel(label);
    }

    public void legend(String position, Font font) {
        CheckExists();
        legend.setItemFont(font);
        legend.setBackgroundPaint(Color.WHITE);
        legend.setFrame(new BlockBorder(Color.BLACK));
        if (position.toLowerCase().equals("northoutside")) {
            legend.setPosition(RectangleEdge.TOP);
            chart.addLegend(legend);
        } else if (position.toLowerCase().equals("eastoutside")) {
            legend.setPosition(RectangleEdge.RIGHT);
            chart.addLegend(legend);
        } else if (position.toLowerCase().equals("southoutside")) {
            legend.setPosition(RectangleEdge.BOTTOM);
            chart.addLegend(legend);
        } else if (position.toLowerCase().equals("westoutside")) {
            legend.setPosition(RectangleEdge.LEFT);
            chart.addLegend(legend);
        } else if (position.toLowerCase().equals("north")) {
            legend.setPosition(RectangleEdge.TOP);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.50,0.98,legend,RectangleAnchor.TOP);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("northeast")) {
            legend.setPosition(RectangleEdge.TOP);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.98,0.98,legend,RectangleAnchor.TOP_RIGHT);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("east")) {
            legend.setPosition(RectangleEdge.RIGHT);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.98,0.50,legend,RectangleAnchor.RIGHT);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("southeast")) {
            legend.setPosition(RectangleEdge.BOTTOM);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.98,0.02,legend,RectangleAnchor.BOTTOM_RIGHT);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("south")) {
            legend.setPosition(RectangleEdge.BOTTOM);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.50,0.02,legend,RectangleAnchor.BOTTOM);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("southwest")) {
            legend.setPosition(RectangleEdge.BOTTOM);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.02,0.02,legend,RectangleAnchor.BOTTOM_LEFT);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("west")) {
            legend.setPosition(RectangleEdge.LEFT);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.02,0.50,legend,RectangleAnchor.LEFT);
            chart.getXYPlot().addAnnotation(ta);
        } else if (position.toLowerCase().equals("northwest")) {
            legend.setPosition(RectangleEdge.TOP);
            XYTitleAnnotation ta = new XYTitleAnnotation(0.02,0.98,legend,RectangleAnchor.TOP_LEFT);
            chart.getXYPlot().addAnnotation(ta);
        }
    }

    public void grid(String xAxis, String yAxis) {
        CheckExists();
        if (xAxis.equalsIgnoreCase("on")){
            chart.getXYPlot().setDomainGridlinesVisible(true);
            chart.getXYPlot().setDomainMinorGridlinesVisible(true);
            chart.getXYPlot().setDomainGridlinePaint(Color.GRAY);
        } else {
            chart.getXYPlot().setDomainGridlinesVisible(false);
            chart.getXYPlot().setDomainMinorGridlinesVisible(false);
        }
        if (yAxis.equalsIgnoreCase("on")){
            chart.getXYPlot().setRangeGridlinesVisible(true);
            chart.getXYPlot().setRangeMinorGridlinesVisible(true);
            chart.getXYPlot().setRangeGridlinePaint(Color.GRAY);
        } else {
            chart.getXYPlot().setRangeGridlinesVisible(false);
            chart.getXYPlot().setRangeMinorGridlinesVisible(false);
        }
    }

    public void saveas(String fileName, int width, int height) {
        CheckExists();
        File file = SemiSim.getUserFile(fileName);
        try {
            ChartUtils.saveChartAsJPEG(file,this.chart,width,height);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void FindColor(String spec, float lineWidth) {
        float dash[] = {5.0f};
        float dot[] = {lineWidth};
        Color color = Color.RED;                    // Default color is red
        Stroke stroke = new BasicStroke(lineWidth); // Default stroke is line
        if (spec.contains("-"))
            stroke = new BasicStroke(lineWidth);
        else if (spec.contains(":"))
            stroke = new BasicStroke(lineWidth, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);
        else if (spec.contains("."))
            stroke = new BasicStroke(lineWidth, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 2.0f, dot, 0.0f);
        if (spec.contains("y"))
            color = Color.YELLOW;
        else if (spec.contains("m"))
            color = Color.MAGENTA;
        else if (spec.contains("c"))
            color = Color.CYAN;
        else if (spec.contains("r"))
            color = Color.RED;
        else if (spec.contains("g"))
            color = Color.GREEN;
        else if (spec.contains("b"))
            color = Color.BLUE;
        else if (spec.contains("k"))
            color = Color.BLACK;
        colors.add(color);
        strokes.add(stroke);
    }
}