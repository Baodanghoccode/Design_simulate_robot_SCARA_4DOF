# SCARA 4-DOF Robot: Trajectory Planning, Dynamic Control & Arduino Implementation

A full-cycle robotics project — from dynamic modeling and control design (MATLAB/Simulink), to mechanical design (SolidWorks), to a Python GUI driving a real Arduino-based prototype.

## Demo

- [Simulation & trajectory visualization](https://youtu.be/2t7z2ixyt3Y)
- [Real hardware motor control](https://youtu.be/OuQhNZBY79Q)

## Pipeline

```
Trajectory planning → Dynamic control (PID) → Hardware (Arduino + motors)
    (.m, cubic poly)      (Simulink, M & C                  ↓
                           matrices)              Visualization (Python GUI)
```

## Repository structure

```
├── Bộ điều khiển robot/
│   ├── Main_Scara4dof.m          # Robot parameters (mass, inertia, link length) + PID gains
│   ├── traj_planning_4dof.m      # Cubic-polynomial joint trajectory generator
│   ├── matrixM.m / matrixC.m     # Mass matrix M(q) and Coriolis matrix C(q, q̇)
│   ├── tinhMa.m / tinhb.m        # Computed-torque control terms
│   ├── RobotScara.slx            # Simulink model: dynamics + PID control loop
│   └── Part*.STEP                # Link geometry (STEP)
├── file vẽ solidwork robot/       # Native SolidWorks CAD (5 parts + assembly)
├── scara_fk_trajectory(...).py    # Python GUI: FK + 3D trajectory simulation (no hardware needed)
├── scara_fk_arduino(...).py       # Python GUI: same visualization + live serial link to Arduino
```

## Technical highlights

**Trajectory planning** — each of the 4 joints follows a smooth 3rd-order polynomial (zero velocity at start/end):

```
q(t) = q0 + (3/T²)(qf - q0)·t² - (2/T³)(qf - q0)·t³
```

**Dynamic control (computed torque / PID)** — uses the robot's own dynamic model to linearize and decouple the system:

```
τ = M(q)·a + C(q, q̇)·q̇ + G(q)
```

with `Kp = diag([5000,5000,5000,5000])`, `Kd = diag([300,300,300,300])`, tuned in Simulink for fast, stable joint tracking.

**Hardware loop** — the Arduino script streams computed joint targets over serial, reads back real position data, and renders the live 3D pose in the same GUI used for pure simulation.

## Tech stack

`MATLAB` · `Simulink` · `SolidWorks` · `Python (Tkinter, Matplotlib, NumPy, PySerial)` · `Arduino`

## Running the simulator

```bash
pip install numpy matplotlib pyserial
python "scara_fk_trajectory(code ve quy dao).py"             # simulation only
python "scara_fk_arduino (code giao dien+ketnoiarduino).py"  # connects to Arduino
```

## Authors

Nguyễn Tiến Bảo · Nguyễn Việt Bách
