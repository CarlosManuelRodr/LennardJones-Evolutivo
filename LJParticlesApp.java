/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.sip.ch08.md2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.controls.XML;
import org.opensourcephysics.display.GUIUtils;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.HistogramFrame;
import org.opensourcephysics.frames.PlotFrame;

class Util {
	public static float Round(float Rval, int Rpl) {
		float p = (float) Math.pow(10, Rpl);
		Rval *= p;
		float tmp = Math.round(Rval);
		return (float) tmp / p;
	}
}

/**
* Comparador de libertad de los sistemas que se usa al ordenar la lista de sistemas 
* en función de su grado de libertad.
*/
class ValorComparator implements Comparator {
	public int compare(Object o1, Object o2) {
		LJParticles a = (LJParticles) o1;
		LJParticles b = (LJParticles) o2;
		return Double.compare(a.getFreedom(), b.getFreedom());
	}

	public boolean equals(Object o) {
		return this == o;
	}
}


/**
 * Crea un conjunto de inviduos en base a los parámetros de un sistema LJParticles inicial.
 */
class Poblacion {
	private java.util.List<LJParticles> individuos;
	
	Poblacion() {
		individuos = new ArrayList<LJParticles>();
	}
	
	Poblacion(long i, LJParticles mdI) {
		individuos = new ArrayList<LJParticles>();

		for (long j = 0; j < i; j++) {
			LJParticles md = new LJParticles();

			md.setDt(mdI.getDt());
			md.setInitialKineticEnergy(mdI.getInitialKineticEnergy());
			md.setLx(mdI.getLx());
			md.setLy(mdI.getLy());
			md.setN(mdI.getN());
			md.setNx(mdI.getNx());
			md.setNy(mdI.getNy());
			md.setRadius(mdI.getRadius());
			md.setRho(mdI.getRho());
			md.setSteps(mdI.getSteps());
			md.setT(mdI.getT());
			md.setAx(mdI.getAx());
			md.setAy(mdI.getAy());
			md.setState(mdI.getState());
			md.setTotalKineticEnergyAccumulator(mdI.getTotalKineticEnergyAccumulator());
			md.setTotalKineticEnergySquaredAccumulator(mdI.getTotalKineticEnergySquaredAccumulator());
			md.setTotalPotentialEnergyAccumulator(mdI.getTotalPotentialEnergyAccumulator());
			md.setVirialAccumulator(mdI.getVirialAccumulator());
			md.getOdeSolver().setStepSize(mdI.getDt());

			individuos.add(md);
		}
	}

	public void evaluar(int steps) {	
		for (LJParticles i : this.getIndividuos()) {
			for (int k = 0; k <= steps; k++) {									
				i.step();					
			}
		}
	}

	public void nuevoIndividuo(Poblacion p) {
		this.individuos.addAll(p.getIndividuos());
	}

	public void nuevoIndividuo(java.util.List<LJParticles> p) {
		this.individuos.addAll(p);
	}

	public java.util.List<LJParticles> getIndividuos() {
		return this.individuos;
	}

	public void setPoblacion(java.util.List<LJParticles> individuos) {
		this.individuos = individuos;
	}
}

/**
 * LJParticlesApp simulates a two-dimensional system of interacting particles
 * via the Lennard-Jones potential.
 * 
 * @author Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.0 revised 03/28/05, 3/29/05
 */
public class LJParticlesApp extends AbstractSimulation {
	LJParticles mdI = new LJParticles();
	LJParticles mdII = new LJParticles();
	PlotFrame temperatureDataI = new PlotFrame("time", "temperature", "Mean temperature system I");
	PlotFrame temperatureDataII = new PlotFrame("time", "temperature", "Mean temperature system II");
	PlotFrame energyDataI = new PlotFrame("time", "Energy", "Mean energy system I");
	PlotFrame energyDataII = new PlotFrame("time", "Energy", "Mean energy system II");
	PlotFrame freedomDataI = new PlotFrame("time", "freedom", "Degree of freedom system I");
	PlotFrame freedomDataII = new PlotFrame("time", "freedom", "Degree of freedom system II");
	DisplayFrame displayI = new DisplayFrame("x", "y", "Lennard-Jones system I");
	DisplayFrame displayII = new DisplayFrame("x", "y", "Lennard-Jones system II");

	Poblacion P;

	private long N = 250;
	private long g = 1250; // 500 generaciones
	private int steps = 100; // 1 segundo
	
	/**
	* @desc Crea una lista con mutaciones de la población. El parámetro que cambia con cada mutación 
	* es la distribución de velocidades.
	* @return Lista con todas las mutaciones.
	*/
	private java.util.List<LJParticles> mutacion() {
		java.util.List<LJParticles> mutaciones = new ArrayList<LJParticles>();
		for (LJParticles i : P.getIndividuos()) {
			LJParticles md = new LJParticles();

			md.setDt(i.getDt());
			md.setInitialKineticEnergy(i.getInitialKineticEnergy());
			md.setInitialConfiguration("rectangular");
			md.setLx(i.getLx());
			md.setLy(i.getLy());
			md.setN(i.getN());
			md.setNx(i.getNx());
			md.setNy(i.getNy());
			md.setRadius(i.getRadius());
			md.setRho(i.getRho());
			md.setSteps(i.getSteps());
			md.setT(i.getT());
			md.setAx(i.getAx());
			md.setAy(i.getAy());
			md.setState(i.getState());
			md.setTotalKineticEnergyAccumulator(i.getTotalKineticEnergyAccumulator());
			md.setTotalKineticEnergySquaredAccumulator(i.getTotalKineticEnergySquaredAccumulator());
			md.setTotalPotentialEnergyAccumulator(i.getTotalPotentialEnergyAccumulator());
			md.setVirialAccumulator(i.getVirialAccumulator());
			md.getOdeSolver().setStepSize(i.getDt());
			md.setVelocitiesEvolution();
			
			mutaciones.add(md);
		}
		
		return mutaciones;
	}
	
	private void buscarSelection(LJParticles mdI) {
		int t = 0;

		P = new Poblacion(N, mdI);

		P.evaluar(steps);
		do {
			
			Poblacion Q = new Poblacion();
			Q.setPoblacion(mutacion());
			
			P.nuevoIndividuo(Q);

			P.evaluar(steps);
			
			java.util.List<LJParticles> resultados = new ArrayList<LJParticles>();
			java.util.List<LJParticles> individuos = new ArrayList<LJParticles>(P.getIndividuos());

			ValorComparator vc = new ValorComparator();
			Collections.sort(individuos, vc);
			long j = 0;

			for (LJParticles i : individuos) {
				if (j < individuos.size() / 2) {
					resultados.add(i);
					j++;
				} else
					break;
			}

			P.setPoblacion(resultados);
			t++;
		} while (t < g);
	}
	
	/**
	 * Initializes the model by reading the number of particles.
	 */
	public void initialize() {
		mdI = new LJParticles();
		mdI.nx = control.getInt("nx"); // number of particles per row
		mdI.ny = control.getInt("ny"); // number of particles per column
		mdI.initialKineticEnergy = control.getDouble("initial kinetic energy per particle");
		mdI.Lx = control.getDouble("Lx");
		mdI.Ly = control.getDouble("Ly");
		mdI.initialConfiguration = control.getString("initial configuration");
		mdI.dt = control.getDouble("dt");

		mdI.initialize();

		for (int k = 0; k < steps; k++) {
			mdI.step();
		}

		displayI.addDrawable(mdI);
		displayI.setPreferredMinMax(0, mdI.Lx, 0, mdI.Ly);

		java.util.List<LJParticles> resultados = new ArrayList<LJParticles>();

		buscarSelection(mdI);
		for (LJParticles i : P.getIndividuos()) {
			resultados.add(i);
		}
		
		mdII = new LJParticles();
		mdII = resultados.get(0);

		displayII.addDrawable(mdII);
		displayII.setPreferredMinMax(0, mdII.Lx, 0, mdII.Ly);
	}

	/**
	 * Does a simulation step and appends data to the views.
	 */
	public void doStep() {

		mdII.step();
		mdII.setVelocitiesEvolution();
		
		energyDataII.append(0, mdII.t, mdII.getMeanEnergy());
		temperatureDataII.append(0, mdII.t, mdII.getMeanTemperature());
		freedomDataII.append(0, mdII.t, mdII.getFreedom());
		
		mdI.step();
		mdI.setVelocitiesEvolution();
		energyDataI.append(0, mdI.t, mdI.getMeanEnergy());
		temperatureDataI.append(0, mdI.t, mdI.getMeanTemperature());
		freedomDataI.append(0, mdI.t, mdI.getFreedom());
	}

	/**
	 * Prints the LJ model's data after the simulation has stopped.
	 */
	public void stop() {
		control.println("Density = " + decimalFormat.format(mdII.rho));
		control.println("Number of time steps = " + mdII.steps);
		control.println("Time step dt = " + decimalFormat.format(mdII.dt));
		control.println("<T>= " + decimalFormat.format(mdII.getMeanTemperature()));
		control.println("<E> = " + decimalFormat.format(mdII.getMeanEnergy()));
		control.println("Heat capacity = " + decimalFormat.format(mdII.getHeatCapacity()));
		control.println("<PA/NkT> = " + decimalFormat.format(mdII.getMeanPressure()));
	}

	/**
	 * Reads adjustable parameters before the program starts running.
	 */
	public void startRunning() {
		mdII.dt = control.getDouble("dt");
		double Lx = control.getDouble("Lx");
		double Ly = control.getDouble("Ly");
		if ((Lx != mdII.Lx) || (Ly != mdII.Ly)) {
			mdII.Lx = Lx;
			mdII.Ly = Ly;
			mdII.computeAcceleration();
			displayII.setPreferredMinMax(0, Lx, 0, Ly);
			resetData();
		}
	}

	/**
	 * Resets the LJ model to its default state.
	 */
	public void reset() {
		control.setValue("nx", 4);
		control.setValue("ny", 4);
		control.setAdjustableValue("Lx", 15.0);
		control.setAdjustableValue("Ly", 10.0);
		control.setValue("initial kinetic energy per particle", 1.0);
		control.setAdjustableValue("dt", 0.01);
		control.setValue("initial configuration", "rectangular");
		enableStepsPerDisplay(true);
		super.setStepsPerDisplay(10); // draw configurations every 10 steps
		displayII.setSquareAspect(true); // so particles will appear as circular disks
		displayI.setSquareAspect(true); // so particles will appear as circular disks
	}

	/**
	 * Resets the LJ model and the data graphs.
	 * 
	 * This method is invoked using a custom button.
	 */
	public void resetData() {
		mdII.resetAverages();
		mdI.resetAverages();
		GUIUtils.clearDrawingFrameData(false); // clears old data from the plot frames
	}

	/**
	 * Returns an XML.ObjectLoader to save and load data for this program.
	 * 
	 * LJParticle data can now be saved using the Save menu item in the control.
	 * 
	 * @return the object loader
	 */
	public static XML.ObjectLoader getLoader() {
		return new LJParticlesLoader();
	}

	/**
	 * Starts the Java application.
	 * 
	 * @param args
	 *            command line parameters
	 */
	public static void main(String[] args) {
		SimulationControl control = SimulationControl.createApp(new LJParticlesApp());
		control.addButton("resetData", "Reset Data");
	}
}

/*
 * Open Source Physics software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.
 * 
 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be
 * released under the GNU GPL license.
 * 
 * This software is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
 * Suite 330, Boston MA 02111-1307 USA or view the license online at
 * http://www.gnu.org/copyleft/gpl.html
 * 
 * Copyright (c) 2007 The Open Source Physics project
 * http://www.opensourcephysics.org
 */
