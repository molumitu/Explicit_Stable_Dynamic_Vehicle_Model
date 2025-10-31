import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from typing import List, Tuple, Dict, Any, Union

# --- 1. Vehicle Model Constants ---
# Vehicle Parameters (Using ALL CAPS for constants)
SPEED = 8.0     # Initial speed (m/s)
RW = 0.325      # Wheel radius (m)
M = 1412.0      # Vehicle mass (m)
A = 1.06        # Distance from CG to front axle (l_f)
B = 1.85        # Distance from CG to rear axle (l_r)
KF = -128915.5  # Front tire cornering stiffness (k_f)
KR = -85943.6   # Rear tire cornering stiffness (k_r)
IZ = 1536.7     # Yaw moment of inertia (I_z)
G = 9.81        # Gravitational acceleration
MU = 0.85       # Road friction coefficient

# --- 2. Plotting Constants and Style Settings ---
PLOT_FONT_SIZE = 18
FONT_FAMILY = 'Times New Roman'

# Apply unified Matplotlib font settings
plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': PLOT_FONT_SIZE,
})

# --- 3. Dynamic Model Functions ---

def dynamic_forward(
    state: np.ndarray, 
    control_input: np.ndarray, 
    dt: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vehicle Dynamic Model using Forward Euler integration (Explicit method).
    
    Args:
        state (x0): Current state vector [X, Y, U, V, Yaw, Omega]
        control_input (u0): Control input vector [a_x, delta]
        dt (T): Time step.

    Returns:
        Next state vector (x1), Lateral force vector [Fy_f, Fy_r], Slip angle vector [alpha_f, alpha_r]
    """
    x0 = state
    u0 = control_input
    T = dt
    
    x1 = np.zeros(len(x0))
    
    # Calculate slip angles (alpha) and lateral forces (Fy) using current state x0
    # Note: abs() is applied as per the original logic, typically used for output/checking.
    alpha1 = abs((x0[3] + A * x0[5] - u0[1] * x0[2]) / x0[2])
    alpha2 = abs((x0[3] - B * x0[5]) / x0[2])
    Fy1 = abs(KF * alpha1)
    Fy2 = abs(KR * alpha2)
    
    # Position: X_next (x1[0])
    x1[0] = x0[0] + T * (x0[2] * np.cos(x0[4]) - x0[3] * np.sin(x0[4]))
    # Position: Y_next (x1[1])
    x1[1] = x0[1] + T * (x0[3] * np.cos(x0[4]) + x0[2] * np.sin(x0[4]))
    # Velocity: U_next (x1[2]) - Longitudinal Velocity
    x1[2] = x0[2] + T * (u0[0] / M / RW + x0[3] * x0[5] - KF * (x0[3] + A * x0[5] - u0[1] * x0[2]) * np.sin(u0[1]) / M / x0[2])
    # Velocity: V_next (x1[3]) - Lateral Velocity
    x1[3] = x0[3] + T * (KF * (x0[3] + A * x0[5] - u0[1] * x0[2]) + KR * (x0[3] - B * x0[5]) - M * x0[2] * x0[2] * x0[5]) / M / x0[2]
    # Yaw Angle: Yaw_next (x1[4])
    x1[4] = x0[4] + T * x0[5]
    # Yaw Rate: Omega_next (x1[5])
    x1[5] = x0[5] + T * (A * KF * (x0[3] + A * x0[5] - u0[1] * x0[2]) - B * KR * (x0[3] - B * x0[5])) / IZ / x0[2]
    
    return x1, np.array([Fy1, Fy2]), np.array([alpha1, alpha2])


def dynamic_ours(
    state: np.ndarray, 
    control_input: np.ndarray, 
    dt: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vehicle Dynamic Model using Backward Euler (Analytic Explicit form).
    
    NOTE: This implementation strictly follows the user's original algebraic formulas for 
    V_next (x1[3]) and Omega_next (x1[5]).
    
    Args:
        state (x0): Current state vector [X, Y, U, V, Yaw, Omega]
        control_input (u0): Control input vector [a_x, delta]
        dt (T): Time step.

    Returns:
        Next state vector (x1), Lateral force vector [Fy_f, Fy_r], Slip angle vector [alpha_f, alpha_r]
    """
    x0 = state
    u0 = control_input
    T = dt
    
    x1 = np.zeros(len(x0))
    
    # Calculate current slip angles and forces (for return value)
    U, V, OMEGA = x0[2], x0[3], x0[5]
    DELTA = u0[1]
    
    if U == 0:
        alpha1 = 0.0
        alpha2 = 0.0
    else:
        alpha1 = abs((V + A * OMEGA - DELTA * U) / U) 
        alpha2 = abs((V - B * OMEGA) / U)
        
    Fy1 = abs(KF * alpha1)
    Fy2 = abs(KR * alpha2)
    
    # --- Explicit States ---
    # Position: X_next (x1[0])
    x1[0] = x0[0] + T * (U * np.cos(x0[4]) - V * np.sin(x0[4]))
    # Position: Y_next (x1[1])
    x1[1] = x0[1] + T * (V * np.cos(x0[4]) + U * np.sin(x0[4]))
    # Velocity: U_next (x1[2])
    x1[2] = x0[2] + T * (u0[0] / M / RW)
    # Yaw Angle: Yaw_next (x1[4])
    x1[4] = x0[4] + T * OMEGA
    
    # --- Implicitly Solved States (User's Original Formulas) ---
    
    # Velocity: V_next (x1[3]) - Lateral Velocity
    NUM_V = (-(A * KF - B * KR) * OMEGA 
             + KF * DELTA * U 
             + M * OMEGA * U**2 
             - M * U * V / T)
             
    DEN_V = KF + KR - M * U / T
    
    if DEN_V == 0:
        x1[3] = V
    else:
        x1[3] = NUM_V / DEN_V

    # Yaw Rate: Omega_next (x1[5])
    NUM_OMEGA = (-IZ * OMEGA * U / T 
                 - (A * KF - B * KR) * V 
                 + A * KF * DELTA * U)
                 
    DEN_OMEGA = ((A**2 * KF + B**2 * KR) - IZ * U / T)
    
    if DEN_OMEGA == 0:
        x1[5] = OMEGA
    else:
        x1[5] = NUM_OMEGA / DEN_OMEGA
    
    return x1, np.array([Fy1, Fy2]), np.array([alpha1, alpha2])


def simulate_trajectory(
    x0: List[float], 
    u_ref: List[List[float]], 
    dt: float, 
    num_steps: int, 
    model_type: str
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Simulates vehicle trajectory using the specified dynamic model.

    Args:
        x0: Initial state.
        u_ref: Sequence of control inputs.
        dt: Time step.
        num_steps: Number of simulation steps.
        model_type: Type of model ('linear_forward' or 'linear_backward').

    Returns:
        State trajectory, Lateral force sequence, Slip angle sequence.
    """
    traj: List[np.ndarray] = [np.array(x0)]
    force: List[np.ndarray] = [np.array([0., 0.])]
    angle: List[np.ndarray] = [np.array([0., 0.])]
    
    if model_type == 'forward':
        dynamic_func = dynamic_forward
    else: # Implies 'ours' based on original logic
        dynamic_func = dynamic_ours
    
    for i in range(num_steps):
        # The result logic follows the original implementation: result = (x1, force, angle)
        result = dynamic_func(traj[i], u_ref[i], dt) 
        traj.append(result[0])
        force.append(result[1])
        angle.append(result[2])
        
    return traj, force, angle


def setup_plot_style(ax: plt.Axes, y_label: str, y_lim: Tuple[float, float]) -> None:
    """Standardizes Matplotlib axis style settings."""
    ax.set_xlabel('Time (s)', fontdict={'family': FONT_FAMILY, 'size': PLOT_FONT_SIZE})
    ax.set_ylabel(y_label, fontdict={'family': FONT_FAMILY, 'size': PLOT_FONT_SIZE})
    ax.set_ylim(y_lim)
    
    # Set tick font properties based on original code's multiple assignments
    for tick in ax.get_xticklabels():
        tick.set_fontname(FONT_FAMILY)
        tick.set_fontsize(PLOT_FONT_SIZE)
    for tick in ax.get_yticklabels():
        tick.set_fontname(FONT_FAMILY)
        tick.set_fontsize(PLOT_FONT_SIZE)

    ax.grid(True)
    bwidth = 1.0
    # Set axis line width
    ax.spines['bottom'].set_linewidth(bwidth)
    ax.spines['left'].set_linewidth(bwidth)
    ax.spines['right'].set_linewidth(bwidth)
    ax.spines['top'].set_linewidth(bwidth)


def run_test(dt: float = 0.1): 
    """
    Executes the vehicle dynamics model simulation test for a given time step (dt).
    """
    T_BENCHMARK = 0.001
    SIMULATION_TIME = 4.0
    NUM_STEPS = int(SIMULATION_TIME / dt)
    
    # Initial state [X, Y, U, V, Yaw, Omega]
    x0: List[float] = [0.0, 0.0, SPEED, 0.0, 0.0, 0.0]

    # --- Control Input: Double step steering ---
    STEER_INPUT = 0.2674 / 2
    step1_len = int(NUM_STEPS / 4)
    # Control input (u_ref) format: [[a_x, delta], ...]
    control_input_sequence = ([[0.0, STEER_INPUT]] * step1_len 
                             + [[0.0, STEER_INPUT * 2]] * (NUM_STEPS - step1_len))
    ur_forw = control_input_sequence
    ur_back = control_input_sequence

    # --- Load Benchmark Data (Ground Truth) ---
    benchmark_file = f"simulink_high_fidelty_groundtruth/simulink_doublestep_v_equals_{int(SPEED)}.mat"
    try:
        benchmark = loadmat(benchmark_file)
        # Flatten and convert to list for consistency with original code's list handling
        w_benchmark: List[float] = benchmark["omega"].flatten().tolist()
        v_benchmark: List[float] = benchmark["v"].flatten().tolist()
        t_benchmark: List[float] = [i * T_BENCHMARK for i in range(len(w_benchmark))]
    except FileNotFoundError:
        print(f"Warning: Benchmark file not found at {benchmark_file}. Skipping plotting benchmark data.")
        w_benchmark, v_benchmark, t_benchmark = [], [], []
    
    BENCHMARK_LABEL = 'Ground truth'
    Y_LIM_V = [-0.05, 1.5]
    Y_LIM_W = [-0.05, 1.5]
        
    # --- Simulation ---
    tra_f, force_f, angle_f = simulate_trajectory(x0, ur_forw, dt, NUM_STEPS, 'forward')
    tra_b, force_b, angle_b = simulate_trajectory(x0, ur_back, dt, NUM_STEPS, 'ours')

    # --- Data Extraction (Using list iteration as in original code) ---
    # Initialize state lists
    v_f, w_f = [], []
    v_b, w_b = [], []
    
    for i in range(NUM_STEPS + 1):
        # Forward Euler data
        v_f.append(tra_f[i][3])
        w_f.append(tra_f[i][5])
        
        # Backward Euler data
        v_b.append(tra_b[i][3])
        w_b.append(tra_b[i][5])
    
    # Time axis for plotted data
    time_points = [i * dt for i in range(NUM_STEPS + 1)]
    
    # --- Plotting ---
    # Create figure and subplots using object-oriented approach (F2.add_subplot)
    fig = plt.figure(figsize=(12, 4), dpi=100)
    ax9 = fig.add_subplot(1, 2, 1) # Lateral velocity (V)
    ax10 = fig.add_subplot(1, 2, 2) # Yaw rate (Omega)
    
    # Subplot 1: Lateral velocity (V)
    ax9.plot(time_points, v_f, '--', color='#FF4500', label='Forward Euler')
    ax9.plot(time_points, v_b, '--', color='#0000CD', label='Our model')
    if t_benchmark: 
        ax9.plot(t_benchmark, v_benchmark, '-', color="#228B22", label=BENCHMARK_LABEL)
    
    setup_plot_style(
        ax9, 
        y_label='Lateral velocity (m/s)', 
        y_lim=Y_LIM_V
    )
    # Specific y-ticks from original code
    ax9.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax9.legend(loc='lower right', prop={'family': FONT_FAMILY, 'size': PLOT_FONT_SIZE})
    
    # Subplot 2: Yaw rate (Omega)
    ax10.plot(time_points, w_f, '--', color='#FF4500', label='Forward Euler')
    ax10.plot(time_points, w_b, '--', color='#0000CD', label='Our model')
    if t_benchmark: 
        ax10.plot(t_benchmark, w_benchmark, '-', color="#228B22", label=BENCHMARK_LABEL)
    
    setup_plot_style(
        ax10, 
        y_label='Yaw rate (rad/s)', 
        y_lim=Y_LIM_W
    )
    # Specific y-ticks from original code
    ax10.set_yticks([0.0, 0.5, 1.0, 1.5])
    
    # Apply tight layout to prevent overlap
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f'results/ts_{int(dt*1000)}ms.pdf', dpi=600, bbox_inches='tight') 
    plt.savefig(f'results/ts_{int(dt*1000)}ms.png', dpi=600, bbox_inches='tight') 
    plt.close(fig) 

# --- 4. Main Execution Logic ---
if __name__ == '__main__':
    # Loop to run tests for different time steps
    for dt_value in [0.1, 0.05, 0.001]:
        run_test(dt_value)