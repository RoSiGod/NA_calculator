import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="Prestress Section Analyzer", layout="wide")

st.title("🏗️ Prestressed RC Section Analyzer")
st.markdown("""
This tool performs a **Serviceability Limit State (SLS)** analysis.
It uses an **Iterative Solver** to find the Neutral Axis (NA) by checking strain compatibility 
and equilibrium for combined Axial Load (P) and Moment (M).
""")

# --- SIDEBAR INPUTS ---
st.sidebar.header("1. Geometry")
b = st.sidebar.number_input("Width (b) [mm]", value=300.0)
h = st.sidebar.number_input("Depth (h) [mm]", value=600.0)

st.sidebar.header("2. Materials")
fc = st.sidebar.number_input("Concrete Strength (fc) [MPa]", value=40.0)
Ec = 4700 * math.sqrt(fc) # ACI approximation
st.sidebar.markdown(f"*Calc. Ec: {Ec:.0f} MPa*")
Es = st.sidebar.number_input("Steel Modulus (Es) [MPa]", value=200000.0)
Eps = st.sidebar.number_input("Strand Modulus (Eps) [MPa]", value=195000.0)

n_s = Es / Ec
n_ps = Eps / Ec

st.sidebar.header("3. Reinforcement")
st.sidebar.subheader("Top Zone (Compression)")
As_top = st.sidebar.number_input("Strand Area (As_top) [mm²]", value=400.0)
d_top = st.sidebar.number_input("Depth to Strands (d_top) [mm]", value=50.0)

st.sidebar.subheader("Bottom Zone (Tension)")
As_bot = st.sidebar.number_input("Rebar Area (As_bot) [mm²]", value=1200.0)
d_bot = st.sidebar.number_input("Depth to Rebar (d_bot) [mm]", value=550.0)

st.sidebar.header("4. Loads")
P = st.sidebar.number_input("Prestress Force (P) [kN]", value=500.0)
M_ext = st.sidebar.number_input("External Moment (M_ext) [kNm]", value=150.0)
load_type = st.sidebar.radio("Moment Type", ["Sagging (Comp. Top)", "Hogging (Tens. Top)"])

# --- SOLVER LOGIC ---

# 1. Standardize Loads (N, Nmm)
P_force = P * 1000 # N (Compression +)
if load_type == "Sagging":
    M_applied = M_ext * 1e6 # Nmm (Compresses Top +)
else:
    M_applied = -M_ext * 1e6 # Nmm (Tensions Top -)

# 2. Check Uncracked State First
# Gross Transformed Properties
Ag = b * h + (n_ps-1)*As_top + (n_s-1)*As_bot
# Moment of Area about Top Fiber
mat_top = (b*h * h/2) + (n_ps-1)*As_top*d_top + (n_s-1)*As_bot*d_bot
y_bar_uncracked = mat_top / Ag

# Inertia about Uncracked Centroid
I_uncracked = (b*h**3)/12 + b*h*(h/2 - y_bar_uncracked)**2 + \
              (n_ps-1)*As_top*(d_top - y_bar_uncracked)**2 + \
              (n_s-1)*As_bot*(d_bot - y_bar_uncracked)**2

# Total Stress Check (Uncracked)
# Shift P and M to the Centroid of the Uncracked Section
# Eccentricity of P relative to Centroid (positive if P is above centroid)
e_p_uncracked = y_bar_uncracked - d_top 
M_total_uncracked = M_applied + (P_force * e_p_uncracked)

# Stress at Bottom Fiber (y is distance from centroid)
y_bot_dist = h - y_bar_uncracked
# Stress = P/A - My/I (Minus because positive M compresses top, so it tensions bottom)
sigma_bot_check = (P_force / Ag) - (M_total_uncracked * y_bot_dist / I_uncracked)

fr = 0.6 * math.sqrt(fc) # Modulus of Rupture
is_cracked = False
if sigma_bot_check < -fr: # Check if tension exceeds limit
    is_cracked = True

# --- FINAL CALCULATION ---

if not is_cracked:
    state = "Uncracked (Elastic)"
    # Calculate stresses using Gross Properties
    sigma_top = (P_force / Ag) + (M_total_uncracked * y_bar_uncracked / I_uncracked)
    sigma_bot = sigma_bot_check
    
    # Calculate NA location (where stress = 0)
    # Using similar triangles between top and bottom stress
    # depth_NA / sigma_top = h / (sigma_top - sigma_bot)
    if sigma_top == sigma_bot:
        NA_final_top = -1 # Pure axial
    else:
        NA_final_top = (sigma_top / (sigma_top - sigma_bot)) * h

else:
    state = "Cracked (Elastic)"
    
    # --- ITERATIVE SOLVER FOR P + M ---
    found_x = False
    best_x = 0
    
    # Iterate through every millimeter of depth to find the true Neutral Axis
    # Logic: Assume NA is at x. Calculate Stress at x. If Stress is 0, we found it.
    for x_try in range(1, int(h)):
        
        # 1. Define Cracked Section Properties for this specific depth x_try
        # Concrete Area (only top part is active)
        Ac = b * x_try
        # Steel Areas (Transformed)
        At_top = (n_ps - 1) * As_top # Compression steel
        At_bot = n_s * As_bot        # Tension steel
        
        A_cr_total = Ac + At_top + At_bot
        
        # Centroid of this specific cracked shape (y_cr) measured from Top
        y_cr = ( (Ac * x_try/2) + (At_top * d_top) + (At_bot * d_bot) ) / A_cr_total
        
        # Inertia of this specific cracked shape about y_cr
        I_cr = (b*x_try**3)/12 + Ac*(x_try/2 - y_cr)**2 + \
               At_top*(d_top - y_cr)**2 + \
               At_bot*(d_bot - y_cr)**2
        
        # 2. Apply Loads to this Centroid
        # Shift P to this new centroid
        e_p_cracked = y_cr - d_top
        M_tot_cracked = M_applied + (P_force * e_p_cracked)
        
        # 3. Calculate Stress at the assumed NA depth (x_try)
        # Distance from centroid to x_try
        dist_centroid_to_assumed_NA = x_try - y_cr
        
        # Stress = P/A + My/I
        # Note: We check stress exactly at x_try. It should be zero.
        sigma_at_boundary = (P_force / A_cr_total) + (M_tot_cracked * dist_centroid_to_assumed_NA / I_cr)
        
        # 4. Convergence Check
        if abs(sigma_at_boundary) < 0.1: # Tolerance of 0.1 MPa
            best_x = x_try
            found_x = True
            
            # Calculate final top/bottom stresses for this valid section
            sigma_top = (P_force / A_cr_total) + (M_tot_cracked * (0 - y_cr) / I_cr)
            sigma_bot = (P_force / A_cr_total) + (M_tot_cracked * (h - y_cr) / I_cr)
            break
            
    if found_x:
        NA_final_top = best_x
    else:
        # HANDLING THE DIVERGENCE
        # If we couldn't find NA inside the section, it usually means 
        # the WHOLE section is in compression (NA > h).
        state = "Fully Compressed (NA > Depth)"
        NA_final_top = h * 1.5 # Just for visualization
        
        # Calculate stress assuming full uncracked section properties again
        # because if NA is outside, the section isn't actually cracked.
        sigma_top = (P_force / Ag) + (M_total_uncracked * y_bar_uncracked / I_uncracked)
        sigma_bot = (P_force / Ag) - (M_total_uncracked * y_bot_dist / I_uncracked)

# --- RESULTS DISPLAY ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Analysis Results")
    if is_cracked:
        st.error(f"**Section State:** {state}")
    else:
        st.success(f"**Section State:** {state}")
    
    st.metric("Neutral Axis Depth (from Top)", f"{NA_final_top:.2f} mm")
    st.metric("Max Concrete Stress (Top)", f"{sigma_top:.2f} MPa")
    st.metric("Stress at Bottom Rebar Level", f"{sigma_bot:.2f} MPa")
    
    if sigma_top > 0.6 * fc:
        st.warning(f"⚠️ High Compressive Stress! ({sigma_top:.1f} > {0.6*fc:.1f} MPa)")

with col2:
    st.subheader("Stress Diagram")
    
    fig, ax = plt.subplots(figsize=(4, 6))
    
    # Draw Concrete Section
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='black', facecolor='#e6e6e6')
    ax.add_patch(rect)
    
    # Draw Strands (Blue)
    circle_ps = patches.Circle((b/2, h - d_top), radius=10, color='blue', label='Prestress')
    ax.add_patch(circle_ps)
    
    # Draw Rebar (Red)
    circle_s = patches.Circle((b/2, h - d_bot), radius=10, color='red', label='Rebar')
    ax.add_patch(circle_s)
    
    # Draw Neutral Axis (Green Dashed)
    na_y = h - NA_final_top
    ax.axhline(y=na_y, color='green', linestyle='--', linewidth=2, label='Neutral Axis')
    
    # Draw Stress Block
    # Base line
    ax.axvline(x=b+40, color='gray', linestyle='-')
    # Top Stress point
    ax.plot([b+40, b+40 + sigma_top*2], [h, h], color='purple')
    # Bot Stress point
    ax.plot([b+40, b+40 + sigma_bot*2], [0, 0], color='purple')
    # Connect
    ax.plot([b+40 + sigma_top*2, b+40 + sigma_bot*2], [h, 0], color='purple', linewidth=2, label='Stress Profile')
    
    ax.set_xlim(-20, b + 150)
    ax.set_ylim(-20, h + 20)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.15))
    ax.set_aspect('equal')
    ax.axis('off')
    st.pyplot(fig)
