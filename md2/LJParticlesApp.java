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
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.frames.HistogramFrame;

/**
* Comparador de libertad de los sistemas que se usa al ordenar la lista de sistemas 
* en función de su grado de libertad.
*/
class ValorComparator implements Comparator<Object> {

	public int compare(Object o1, Object o2) {
		Individuo a = (Individuo) o1;
		Individuo b = (Individuo) o2;

		return Double.compare(a.getFreedom(), b.getFreedom());
	}

	public boolean equals(Object o) {
		return this == o;
	}
}

/**
 * Contiene sistema LJ junto con su valor de libertad y cambio de
 * energía del sistema con evolución guiada y el sistema sin evolución
 * guiada (mdI).
 */
class Individuo {
	private LJParticles lj;
	private double value;
	private double valueMdI;
	private double freedom;
	private double deltaF; // Cambio de libertad.

	public double getFreedom() {
		return freedom;
	}

	public void setFreedom(double freedom) {
		this.freedom = freedom;
	}
	
	public double getDeltaFreedom() {
		return deltaF;
	}

	public void setDeltaFreedom(double deltaF) {
		this.deltaF = deltaF;
	}

	public double getValueMdI() {
		return valueMdI;
	}

	public void setValueMdI(double valueMdI) {
		this.valueMdI = valueMdI;
	}
	
	public void setValue(double value) {
		this.value = value;
	}

	Individuo(LJParticles lj) {
		this.lj = lj;
		value = 0d;
		freedom = 0d;
	}

	public double getValue() {
		return value;
	}

	public LJParticles getLj() {
		return lj;
	}

	public void setLj(LJParticles lj) {
		this.lj = lj;
	}
}

/**
* Crea un conjunto de inviduos en base a los parámetros de un sistema LJParticles inicial.
*/
class Poblacion {
	private java.util.List<Individuo> individuos;
	
	Poblacion() {
		individuos = new ArrayList<Individuo>();
	}
	
	Poblacion(long i, LJParticles mdI) {
		individuos = new ArrayList<Individuo>();

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

			Individuo I = new Individuo(md);
			I.setFreedom(mdI.getMeanFreedom());
			I.setValue(mdI.getMeanEnergy());
			
			individuos.add(I);
		}
	}

	/**
	 * Realiza la evolución temporal de la población y del sistema mdI guardando
	 * los cambio de energía.
	 * 
	 * @param steps Número de pasos.
	 * @param mdI Sistema LJ.
	 */
	public void evaluar(int steps, LJParticles mdI) {
		
		double Ei = mdI.getEnergy();
		
		for(int k = 0; k <= steps; k++)
		{									
			mdI.step();				
		}
		
		double Ej = mdI.getEnergy();
		double mdIBQ = (-1)*(Ei - Ej);
		double Fi, Fj;
		
		for (Individuo i : this.getIndividuos()) {
			
			Ei = i.getLj().getEnergy();
			Fi = i.getLj().getFreedom();
			
			for(int k = 0; k <= steps; k++)
			{									
				i.getLj().step();					
			}
			
			Ej = i.getLj().getEnergy();
			Fj = i.getLj().getFreedom();
			
			double BQ = (-1)*(Ei - Ej); 
				
			i.setValue(i.getValue() + BQ);
			i.setFreedom(i.getLj().getMeanFreedom());
			i.setDeltaFreedom(Fj - Fi);
			i.setValueMdI(i.getValueMdI() + mdIBQ);
		}
	}

	public void nuevoIndividuo(Poblacion p) {
		this.individuos.addAll(p.getIndividuos());
	}

	public void nuevoIndividuo(java.util.List<Individuo> p) {
		this.individuos.addAll(p);
	}

	public java.util.List<Individuo> getIndividuos() {
		return this.individuos;
	}

	public void setPoblacion(java.util.List<Individuo> individuos) {
		this.individuos = individuos;
	}
}

/**
 * LJParticlesApp simula un sistema bidimensional de partículas que interactúan
 * entre ellas a través de un potencial de Lennard-Jones y se utiliza un algoritmo
 * evolutivo para seleccionar las mutaciones que presentan menor grado de libertad.
 * 
 * @author Jan Tobochnik, Wolfgang Christian, Harvey Gould, Samuel Galván
 * @version 1.0 revised 03/28/05, 3/29/05
 */
public class LJParticlesApp extends AbstractSimulation {
	LJParticles mdI;
	LJParticles mdII;
	LJParticles mdIII;
	PlotFrame temperatureDataI;
	PlotFrame temperatureDataII;
	PlotFrame energyDataI;
	PlotFrame energyDataII;
	PlotFrame freedomDataI;
	PlotFrame freedomDataII;
	PlotFrame meanDeltaFreedomIII;
	DisplayFrame displayI;
	DisplayFrame displayII;
	DisplayFrame displayIII;
	PlotFrame freedomTemperatureData;
	HistogramFrame deltaEDist;
	HistogramFrame deltaFreedDist;

	Poblacion P, PmdIII;
	private long N;
	private long g;
	private int steps;
	private boolean simular_sin_mutar, simular_mutado, simular_mutaciones, ver_dist_cambio_energ;
	
	/**
	* Crea una lista con mutaciones de la población. El parámetro que cambia con cada mutación 
	* es la distribución de velocidades.
	* @return Lista con todas las mutaciones.
	*/
	private java.util.List<Individuo> mutacion(Poblacion pmut){
		java.util.List<Individuo> mutaciones = new ArrayList<Individuo>();
		for (Individuo i : pmut.getIndividuos()) {
			LJParticles md = new LJParticles();

			md.setDt(i.getLj().getDt());
			md.setInitialKineticEnergy(i.getLj().getInitialKineticEnergy());
			md.setInitialConfiguration("rectangular");
			md.setLx(i.getLj().getLx());
			md.setLy(i.getLj().getLy());
			md.setN(i.getLj().getN());
			md.setNx(i.getLj().getNx());
			md.setNy(i.getLj().getNy());
			md.setRadius(i.getLj().getRadius());
			md.setRho(i.getLj().getRho());
			md.setSteps(i.getLj().getSteps());
			md.setT(i.getLj().getT());
			md.setAx(i.getLj().getAx());
			md.setAy(i.getLj().getAy());
			md.setState(i.getLj().getState());
			md.setTotalKineticEnergyAccumulator(i.getLj().getTotalKineticEnergyAccumulator());
			md.setTotalKineticEnergySquaredAccumulator(i.getLj().getTotalKineticEnergySquaredAccumulator());
			md.setTotalPotentialEnergyAccumulator(i.getLj().getTotalPotentialEnergyAccumulator());
			md.setVirialAccumulator(i.getLj().getVirialAccumulator());
			md.getOdeSolver().setStepSize(i.getLj().getDt());
			
			md.setVelocitiesEvolution();
			
			Individuo I = new Individuo(md);			
			I.setFreedom(i.getFreedom());
			I.setValue(i.getValue());
			I.setValueMdI(i.getValueMdI());
			
			mutaciones.add(I);
		}
		
		return mutaciones;
	}
	
	/**
	* Evoluciona MDI buscando las mutaciones que tienen menor grado de libertad.
	**/
	private void buscarSelection(LJParticles mdI) {
		int t = 0;
		P = new Poblacion(N, mdI);
		P.evaluar(steps, mdI);
		
		do {
			Poblacion Q = new Poblacion();
			Q.setPoblacion(mutacion(P));
			P.nuevoIndividuo(Q);
			P.evaluar(steps, mdI);
			
			java.util.List<Individuo> resultados = new ArrayList<Individuo>();
			java.util.List<Individuo> individuos = new ArrayList<Individuo>(P.getIndividuos());

			ValorComparator vc = new ValorComparator();
			Collections.sort(individuos, vc);
			long j = 0;

			for (Individuo i : individuos) {
				if (j < individuos.size() / 2) {
					resultados.add(i);
					j++;
				} else
					break;
			}

			P.setPoblacion(resultados);
			
			// Seleccionamos el individuo con menos complejidad y graficamos su consumo energético.
			Individuo fit = P.getIndividuos().get(0);
			freedomTemperatureData.append(0, fit.getFreedom(), fit.getValueMdI() -  fit.getValue());
			t++;
		} while (t < g);
	}
	
	/**
	 * Ejecuta el algoritmo evolutivo.
	 **/
	public void initialize() {	
		// Copia los parámetros de la interfaz gráfica.
		g = control.getInt("Generaciones");
		steps = control.getInt("Pasos por generación");
		N = control.getInt("Tamaño de la población");
		simular_sin_mutar = control.getBoolean("Simular sistema sin mutar");
		simular_mutado = control.getBoolean("Simular sistema mutado");
		simular_mutaciones = control.getBoolean("Simular mutaciones");
		ver_dist_cambio_energ = control.getBoolean("Ver distribución de cambio energético");
		
		if (simular_sin_mutar || simular_mutado)
		{
				super.setStepsPerDisplay(10);  // Dibuja la configuración cada 10 pasos.
				mdI = new LJParticles();
				mdI.nx = control.getInt("nx"); // Número de partículas por fila.
				mdI.ny = control.getInt("ny"); // Número de partículas por columna.
				mdI.initialKineticEnergy = control.getDouble("Energía cinética inicial por partícula");
				mdI.Lx = control.getDouble("Lx");
				mdI.Ly = control.getDouble("Ly");
				mdI.initialConfiguration = control.getString("Configuración inicial");
				mdI.dt = control.getDouble("dt");
		
				// Crea sistema LJ mdI y lo evoluciona un número 'step' de pasos.
				mdI.initialize();
		
				for (int k = 0; k < steps; k++) {
					mdI.step();
				}
		
				if (simular_sin_mutar)
				{
					displayI = new DisplayFrame("x", "y", "Sistema Lennard-Jones I");
					displayI.setSquareAspect(true);
					temperatureDataI = new PlotFrame("Tiempo", "Temperatura", "Temperatura promedio en sistema I");
					energyDataI = new PlotFrame("Tiempo", "Energía", "Energía promedio en sistema I");
					freedomDataI = new PlotFrame("Tiempo", "Libertad", "Grado de libertad en sistema I");
					
					displayI.addDrawable(mdI);
					displayI.setPreferredMinMax(0, mdI.Lx, 0, mdI.Ly);
				}
		
				if (simular_mutado)
				{
					displayII = new DisplayFrame("x", "y", "Sistema Lennard-Jones II");
					displayII.setSquareAspect(true);
					temperatureDataII = new PlotFrame("Tiempo", "Temperatura", "Temperatura promedio en sistema II");
					energyDataII = new PlotFrame("Tiempo", "Energía", "Energía promedio en sistema II");
					freedomDataII = new PlotFrame("Tiempo", "Libertad", "Grado de libertad en sistema II");
					freedomTemperatureData = new PlotFrame("Libertad", "Calor disipado",
							"Relación entre calor disipado vs libertad");
					freedomTemperatureData.setLogScale(true, true);
					
					java.util.List<Individuo> resultados = new ArrayList<Individuo>();
			
					// Evoluciona mdI seleccionando las mutaciones con menor grado de libertad.
					buscarSelection(mdI);
					for (Individuo i : P.getIndividuos()) {
						resultados.add(i);
						
					// Crea el sistema mdII a partir de la mutación con menor libertad de mdI.
					mdII = new LJParticles();
					mdII = resultados.get(0).getLj();
					displayII.addDrawable(mdII);
					displayII.setPreferredMinMax(0, mdII.Lx, 0, mdII.Ly);
				}
			}
		}
		
		if (simular_mutaciones)
		{
			super.setStepsPerDisplay(1);
			displayIII = new DisplayFrame("x", "y", "Sistema Lennard-Jones III");
			displayIII.setSquareAspect(true);
			deltaFreedDist = new HistogramFrame("Cambio de libertad", "Frecuencia",
					"Distribución de cambios de libertad en el sistema III");
			meanDeltaFreedomIII = new PlotFrame("Tiempo", "Cambio de libertad promedio",
					"Cambio de libertad promedio en sistema III");
			
			// El sistema mdIII se usa para calcular la distribución de los valores de libertad.
			mdIII = new LJParticles();
			mdIII.nx = control.getInt("nx");
			mdIII.ny = control.getInt("ny");
			mdIII.initialKineticEnergy = control.getDouble("Energía cinética inicial por partícula");
			mdIII.Lx = control.getDouble("Lx");
			mdIII.Ly = control.getDouble("Ly");
			mdIII.initialConfiguration = control.getString("Configuración inicial");
			mdIII.dt = control.getDouble("dt");
			mdIII.initialize();
			PmdIII = new Poblacion(N, mdIII);
			PmdIII.evaluar(steps, mdIII);
			displayIII.addDrawable(mdIII);
			displayIII.setPreferredMinMax(0, mdIII.Lx, 0, mdIII.Ly);
			deltaFreedDist.setBinWidth(4.0);
		}
		
		if (ver_dist_cambio_energ)
		{
			// El sistema mdIV se usa para calcular los valores de la distribución de cambio
			// de energía.
			deltaEDist = new HistogramFrame("Cambio de energía", "Frecuencia",
					"Distribución de cambios de energía");
			LJParticles mdIV = new LJParticles();
			mdIV.nx = control.getInt("nx");
			mdIV.ny = control.getInt("ny");
			mdIV.initialKineticEnergy = control.getDouble("Energía cinética inicial por partícula");
			mdIV.Lx = control.getDouble("Lx");
			mdIV.Ly = control.getDouble("Ly");
			mdIV.initialConfiguration = control.getString("Configuración inicial");
			mdIV.dt = control.getDouble("dt");
			mdIV.initialize();
			deltaEDist.setBinWidth(4.0);
			deltaFreedDist.setBinWidth(4.0);
			
			for (int i = 0; i < g; i++)
			{
				double Ei = mdIV.getEnergy();
				for (int j = 0; j <= steps; j++)
				{									
					mdIV.step();	
					mdIV.setVelocitiesEvolution();
				}
				double Ej = mdIV.getEnergy();
				
				deltaEDist.append((-1)*(Ei - Ej));
			}
		}
	}

	/**
	 * Efectúa un paso de la simulación y actualiza las gráficas.
	 */
	public void doStep() {
		if (simular_sin_mutar)
		{
			// Evoluciona mdI y actualiza gráficas.
			mdI.step();
			mdI.setVelocitiesEvolution();
			energyDataI.append(0, mdI.t, mdI.getMeanEnergy());
			temperatureDataI.append(0, mdI.t, mdI.getMeanTemperature());
			freedomDataI.append(0, mdI.t, mdI.getFreedom());
		}
		
		if (simular_mutado)
		{
			// Evoluciona mdII y actualiza gráficas.
			mdII.step();
			mdII.setVelocitiesEvolution();
			energyDataII.append(0, mdII.t, mdII.getMeanEnergy());
			temperatureDataII.append(0, mdII.t, mdII.getMeanTemperature());
			freedomDataII.append(0, mdII.t, mdII.getFreedom());
		}
		
		if (simular_mutaciones)
		{
			// Evoluciona mdIII
			Poblacion Q = new Poblacion();
			Q.setPoblacion(mutacion(PmdIII));
			PmdIII.nuevoIndividuo(Q);
			PmdIII.evaluar(steps, mdIII);
			java.util.List<Individuo> resultados = new ArrayList<Individuo>();
			java.util.List<Individuo> individuos = new ArrayList<Individuo>(PmdIII.getIndividuos());
			
			ValorComparator vc = new ValorComparator();
			Collections.sort(individuos, vc);
			long j = 0;
	
			deltaFreedDist.setBinWidth(15.0);
			
			double meanFreedom = 0;
			for (Individuo i : individuos) {
				if (j < individuos.size() / 2) {
					resultados.add(i);
					deltaFreedDist.append(i.getDeltaFreedom());
					meanFreedom += i.getDeltaFreedom();
					j++;
				} else
					break;
			}
			meanFreedom /= individuos.size();
			meanDeltaFreedomIII.append(0, mdIII.t, meanFreedom);
			
			PmdIII.setPoblacion(resultados);
			displayIII.clearDrawables();
			displayIII.addDrawable(PmdIII.getIndividuos().get(0).getLj());
		}
	}

	/**
	 * Escribe la información del sistema LJ al terminar la simulación.
	 */
	public void stop() {
		if (simular_mutado)
		{
			control.println("Density = " + decimalFormat.format(mdII.rho));
			control.println("Number of time steps = " + mdII.steps);
			control.println("Time step dt = " + decimalFormat.format(mdII.dt));
			control.println("<T>= " + decimalFormat.format(mdII.getMeanTemperature()));
			control.println("<E> = " + decimalFormat.format(mdII.getMeanEnergy()));
			control.println("Heat capacity = " + decimalFormat.format(mdII.getHeatCapacity()));
			control.println("<PA/NkT> = " + decimalFormat.format(mdII.getMeanPressure()));
		}
	}

	/**
	 * Copia los parámetros del control antes de comenzar la ejecución.
	 */
	public void startRunning() {
		if (simular_mutado)
		{
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
	}

	/**
	 * Resetea el modelo de LJ a sus valores por defecto.
	 */
	public void reset() {
		control.setValue("nx", 4);
		control.setValue("ny", 4);
		control.setValue("Generaciones", 1250);
		control.setValue("Pasos por generación", 100);
		control.setValue("Tamaño de la población", 250);
		control.setAdjustableValue("Lx", 15.0);
		control.setAdjustableValue("Ly", 10.0);
		control.setValue("Energía cinética inicial por partícula", 1.0);
		control.setAdjustableValue("dt", 0.01);
		control.setValue("Configuración inicial", "rectangular");
		control.setValue("Simular sistema sin mutar", true);
		control.setValue("Simular sistema mutado", true);
		control.setValue("Simular mutaciones", false);
		control.setValue("Ver distribución de cambio energético", false);
		enableStepsPerDisplay(true);
	}

	/**
	 * Resetea el modelo de LJ y las gráficas.
	 */
	public void resetData() {
		if (simular_mutado)
			mdII.resetAverages();
		
		if (simular_sin_mutar)
			mdI.resetAverages();
		
		GUIUtils.clearDrawingFrameData(false); // clears old data from the plot
											   // frames
	}

	/**
	 * Devuelve el objeto XML.ObjectLoader para cargar y guardar datos.
	 */
	public static XML.ObjectLoader getLoader() {
		return new LJParticlesLoader();
	}

	/**
	 * Comienza la ejecución del programa.
	 * 
	 * @param args
	 *            Argumentos de la línea de comandos.
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
