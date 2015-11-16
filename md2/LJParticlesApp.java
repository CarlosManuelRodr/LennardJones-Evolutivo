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


/**
* Comparador de libertad de los sistemas que se usa al ordenar la lista de sistemas 
* en funci�n de su grado de libertad.
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
 * Contiene sistema LJ junto con su valor de libertad y cambio de energ�a.
 */
class Individuo {
	private LJParticles lj;
	private double value;
	private double valueMdI;
	private double freedom;

	public double getFreedom() {
		return freedom;
	}

	public void setFreedom(double freedom) {
		this.freedom = freedom;
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
* Crea un conjunto de inviduos en base a los par�metros de un sistema LJParticles inicial.
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
			md.setTotalKineticEnergyAccumulator(mdI
					.getTotalKineticEnergyAccumulator());
			md.setTotalKineticEnergySquaredAccumulator(mdI
					.getTotalKineticEnergySquaredAccumulator());
			md.setTotalPotentialEnergyAccumulator(mdI
					.getTotalPotentialEnergyAccumulator());
			md.setVirialAccumulator(mdI.getVirialAccumulator());
			md.getOdeSolver().setStepSize(mdI.getDt());

			Individuo I = new Individuo(md);
			I.setFreedom(mdI.getMeanFreedom());
			I.setValue(mdI.getMeanEnergy());
			
			individuos.add(I);
		}
	}

	/**
	 * Realiza la evoluci�n temporal de la poblaci�n y del sistema mdI guardando
	 * los cambio de energ�a.
	 * 
	 * @param steps N�mero de pasos.
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
		
		for (Individuo i : this.getIndividuos()) {
			
			Ei = i.getLj().getEnergy();
			
			for(int k = 0; k <= steps; k++)
			{									
				i.getLj().step();					
			}
			
			Ej = i.getLj().getEnergy();
			
			double BQ = (-1)*(Ei - Ej); 
				
			i.setValue(i.getValue() + BQ);
			i.setFreedom(i.getLj().getMeanFreedom());
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
 * LJParticlesApp simula un sistema bidimensional de part�culas que interact�an
 * entre ellas a trav�s de un potencial de Lennard-Jones y se utiliza un algoritmo
 * evolutivo para seleccionar las mutaciones que presentan menor grado de libertad.
 * 
 * @author Jan Tobochnik, Wolfgang Christian, Harvey Gould, Samuel Galv�n
 * @version 1.0 revised 03/28/05, 3/29/05
 */
public class LJParticlesApp extends AbstractSimulation {
	LJParticles mdI = new LJParticles();
	LJParticles mdII = new LJParticles();
	PlotFrame temperatureDataI = new PlotFrame("Tiempo", "Temperatura", "Temperatura promedio en sistema I");
	PlotFrame temperatureDataII = new PlotFrame("Tiempo", "Temperatura", "Temperatura promedio en sistema II");
	PlotFrame energyDataI = new PlotFrame("Tiempo", "Energ�a", "Energ�a promedio en sistema I");
	PlotFrame energyDataII = new PlotFrame("Tiempo", "Energ�a", "Energ�a promedio en sistema II");
	PlotFrame freedomDataI = new PlotFrame("Tiempo", "Libertad", "Grado de libertad en sistema I");
	PlotFrame freedomDataII = new PlotFrame("Tiempo", "Libertad", "Grado de libertad en sistema II");
	DisplayFrame displayI = new DisplayFrame("x", "y", "Sistema Lennard-Jones I");
	DisplayFrame displayII = new DisplayFrame("x", "y", "Sistema Lennard-Jones II");
	PlotFrame freedomTemperatureData = new PlotFrame("Libertad", "Calor promedio",
			"Relaci�n entre calor promedio vs libertad");

	Poblacion P;
	private long N;
	private long g;
	private int steps;
	
	/**
	* Crea una lista con mutaciones de la poblaci�n. El par�metro que cambia con cada mutaci�n 
	* es la distribuci�n de velocidades.
	* @return Lista con todas las mutaciones.
	*/
	private java.util.List<Individuo> mutacion() {
		java.util.List<Individuo> mutaciones = new ArrayList<Individuo>();
		for (Individuo i : P.getIndividuos()) {
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
			
			md.setVelocitiesEvolution(); // Cambia la velocidad de las part�culas con 10% de probabilidad.
			
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
			// Realiza mutaciones sobre la poblaci�n P, las a�ade a la lista de poblaci�n y evoluciona.
			Poblacion Q = new Poblacion();
			Q.setPoblacion(mutacion());
			
			P.nuevoIndividuo(Q);
			P.evaluar(steps, mdI);
			
			java.util.List<Individuo> resultados = new ArrayList<Individuo>();
			java.util.List<Individuo> individuos = new ArrayList<Individuo>(P.getIndividuos());

			// Ordena la lista de individuos seg�n su valor de libertad.
			ValorComparator vc = new ValorComparator();
			Collections.sort(individuos, vc);
			
			// Descarta los peores resultados.
			long j = 0;
			for (Individuo i : individuos) {
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
	 * Ejecuta el algoritmo evolutivo.
	 **/
	public void initialize() {
		// Copia los par�metros de la interfaz gr�fica.
		g = control.getInt("Generaciones");
		steps = control.getInt("Pasos por generaci�n");
		N = control.getInt("Tama�o de la poblaci�n");
		
		mdI = new LJParticles();
		mdI.nx = control.getInt("nx"); // N�mero de part�culas por fila.
		mdI.ny = control.getInt("ny"); // N�mero de part�culas por columna.
		mdI.initialKineticEnergy = control.getDouble("Energ�a cin�tica inicial por part�cula");
		mdI.Lx = control.getDouble("Lx");
		mdI.Ly = control.getDouble("Ly");
		mdI.initialConfiguration = control.getString("Configuraci�n inicial");
		mdI.dt = control.getDouble("dt");

		// Crea sistema LJ mdI y lo evoluciona un n�mero 'step' de pasos.
		mdI.initialize();

		for (int k = 0; k < steps; k++) {
			mdI.step();
		}

		displayI.addDrawable(mdI);
		displayI.setPreferredMinMax(0, mdI.Lx, 0, mdI.Ly);

		java.util.List<Individuo> resultados = new ArrayList<Individuo>();

		// Evoluciona mdI seleccionando las mutaciones con menor grado de libertad.
		buscarSelection(mdI);
		for (Individuo i : P.getIndividuos()) {
			resultados.add(i);
		}

		for (Individuo i : resultados) {
			freedomTemperatureData.append(0, i.getFreedom(), i.getValue());
		}
		
		// Crea el sistema mdII a partir de la mutaci�n con menor libertad de mdI.
		mdII = new LJParticles();
		mdII = resultados.get(0).getLj();

		displayII.addDrawable(mdII);
		displayII.setPreferredMinMax(0, mdII.Lx, 0, mdII.Ly);
	}

	/**
	 * Efect�a un paso de la simulaci�n y actualiza las gr�ficas.
	 */
	public void doStep() {
		// Evoluciona mdI y actualiza gr�ficas.
		mdI.step();
		mdI.setVelocitiesEvolution();
		energyDataI.append(0, mdI.t, mdI.getMeanEnergy());
		temperatureDataI.append(0, mdI.t, mdI.getMeanTemperature());
		freedomDataI.append(0, mdI.t, mdI.getFreedom());
		
		// Evoluciona mdII y actualiza gr�ficas.
		mdII.step();
		mdII.setVelocitiesEvolution();
		energyDataII.append(0, mdII.t, mdII.getMeanEnergy());
		temperatureDataII.append(0, mdII.t, mdII.getMeanTemperature());
		freedomDataII.append(0, mdII.t, mdII.getFreedom());
	}

	/**
	 * Escribe la informaci�n del sistema LJ al terminar la simulaci�n.
	 */
	public void stop() {
		control.println("Densidad = " + decimalFormat.format(mdII.rho));
		control.println("N�mero de pasos = " + mdII.steps);
		control.println("Intervalo temporal de paso = " + decimalFormat.format(mdII.dt));
		control.println("<T>= " + decimalFormat.format(mdII.getMeanTemperature()));
		control.println("<E> = " + decimalFormat.format(mdII.getMeanEnergy()));
		control.println("Capacidad calor�fica = " + decimalFormat.format(mdII.getHeatCapacity()));
		control.println("<PA/NkT> = " + decimalFormat.format(mdII.getMeanPressure()));
	}

	/**
	 * Copia los par�metros del control antes de comenzar la ejecuci�n.
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
	 * Resetea el modelo de LJ a sus valores por defecto.
	 */
	public void reset() {
		control.setValue("nx", 4);
		control.setValue("ny", 4);
		control.setValue("Generaciones", 1250);
		control.setValue("Pasos por generaci�n", 100);
		control.setValue("Tama�o de la poblaci�n", 250);
		control.setAdjustableValue("Lx", 15.0);
		control.setAdjustableValue("Ly", 10.0);
		control.setValue("Energ�a cin�tica inicial por part�cula", 1.0);
		control.setAdjustableValue("dt", 0.01);
		control.setValue("Configuraci�n inicial", "rectangular");
		enableStepsPerDisplay(true);
		super.setStepsPerDisplay(10); // draw configurations every 10 steps
		displayII.setSquareAspect(true); // so particles will appear as circular disks
		displayI.setSquareAspect(true); // so particles will appear as circular disks
	}

	/**
	 * Resetea el modelo de LJ y las gr�ficas.
	 */
	public void resetData() {
		mdII.resetAverages();
		mdI.resetAverages();
		GUIUtils.clearDrawingFrameData(false); // clears old data from the plot frames
	}

	/**
	 * Devuelve un XML.ObjectLoader para guardar y cargar los datos de este programa.
	 * Los datos de LJParticle se pueden guardar usando el menu.
	 * 
	 * @return El ObjectLoader.
	 */
	public static XML.ObjectLoader getLoader() {
		return new LJParticlesLoader();
	}

	/**
	 * Comienza la ejecuci�n de la aplicaci�n.
	 * @param args argumentos de la l�nea de comandos.
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
