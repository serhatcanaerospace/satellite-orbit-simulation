# 🌍 Satellite Orbit Simulation and Analysis Tool

A MATLAB-based tool designed to simulate satellite orbits, analyze perturbation effects, and predict collision risks using real-world TLE data. This project is perfect for students, researchers, and professionals in aerospace engineering.

---

## 🚀 Key Features

- **📊 TLE Data Integration:** Processes Two-Line Element (TLE) data to extract orbital parameters.
- **🌌 3D Orbit Visualization:** Simulates and visualizes satellite trajectories in 3D.
- **📈 Perturbation Analysis:** Models long-term effects of J2 perturbation, atmospheric drag, and solar radiation pressure.
- **⚠️ Collision Risk Prediction:** Calculates proximity between satellites and predicts collision risks.
- **🛰️ Multi-Satellite Simulation:** Simulates constellations like Starlink for advanced analysis.

---

## 🛠️ How It Works

1. **Load TLE Data:** Start with real TLE data, such as the provided `last-30-days.txt`.
2. **Process Orbital Parameters:** Extract inclination, eccentricity, RAAN, and more.
3. **Visualize Orbits:** Generate 3D plots of satellite orbits.
4. **Analyze Perturbations:** Observe long-term changes in orbital parameters.
5. **Predict Risks:** Identify potential satellite collisions.

---

## 🖥️ Getting Started

### Requirements
- **MATLAB:** Version R2021a or later
- **TLE Data File:** Example provided in the repository (`last-30-days.txt`)

### Clone the Repository
```bash
git clone https://github.com/seroserosero123123/satellite-orbit-simulation.git
Run in MATLAB
Open main_script.m in MATLAB.
Ensure the TLE data file (last-30-days.txt) is in the same directory.
Run the script to process data and generate visualizations.
📂 Project Structure
main_script.m: Main script to execute the project.
readTLE.m: Reads and processes TLE data.
processTLE.m: Converts TLE data into orbital parameters.
last-30-days.txt: Example TLE data for analysis.
outputs/: Folder containing visualization results.
🌌 Visualization Examples
3D Orbit Visualization


Perturbation Analysis


👨‍💻 About the Author
Serhat Dalmış
Aerospace Engineering Student passionate about orbital mechanics, satellite systems, and computational tools for space exploration.

🌐 LinkedIn
📧 serhatcandalmis@gmail.com
📊 GitHub Stats

🤝 Contributions
Contributions, issues, and feature requests are welcome! Feel free to fork this repository or open an issue to suggest improvements.

🛰️ Future Plans
Expand Perturbation Modeling: Add more detailed atmospheric drag and solar radiation effects.
Interactive User Interface: Build a MATLAB GUI for easier use.
Exportable Reports: Automate generation of PDF reports from simulation results.
