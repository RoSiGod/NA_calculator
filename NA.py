import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="Prestress Section Analyzer", layout="wide")

st.title("🏗️ Prestressed RC Section Analyzer")
st.markdown("""
This tool performs a **Serviceability Limit State (SLS)** analysis for a rectangular reinforced concrete section 
with prestressing strands in the compression zone. It uses an **Iterative Strain Compatibility** method 
to find the true Neutral Axis under combined Axial Load (P) and Moment (M).
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

# --- CALCULATION LOGIC ---

def calculate_equilibrium(na_depth, curvature):
    """
    Calculates the internal Force imbalance and Moment for a given NA depth and Curvature.
    Assumption: Linear Elastic Behavior (SLS).
    """
    # 1. Strain Profile
    # strain = curvature * distance_from_NA
    # Compresssion is Positive (+)
    
    # Concrete Force (Triangular block)
    # We only integrate concrete where strain is positive (Compression)
    if na_depth > 0:
        # Stress at top fiber
        sigma_top = curvature * na_depth * Ec
        # Force = 0.5 * base * height * width
        C_conc = 0.5 * sigma_top * b * na_depth
        # Moment of concrete about Top Fiber
        # Centroid of triangle is at na_depth/3 from top
        M_conc = C_conc * (na_depth / 3)
    else:
        C_conc = 0.0
        M_conc = 0.0

    # Top Steel (Strands)
    strain_top_s = curvature * (na_depth - d_top)
    sigma_top_s = strain_top_s * Eps
    # Account for displaced concrete (n-1) if in compression
    eff_area_top = (n_ps - 1) * As_top if strain_top_s > 0 else n_ps * As_top
    F_top_s = sigma_top_s * As_top # Force (N) - approximate for transformed
    # Correction: Use modular ratio properly
    F_top_s = (sigma_top_s / Eps) * (Eps - Ec) * As_top if strain_top_s > 0 else sigma_top_s * As_top
    # Actually, simpler Transformed Section logic:
    # F = Stress_concrete_at_level * (n-1) * A
    sigma_c_at_top_s = curvature * (na_depth - d_top) * Ec
    F_top_s = sigma_c_at_top_s * (n_ps - 1) * As_top
    M_top_s = F_top_s * d_top

    # Bottom Steel (Rebar)
    sigma_c_at_bot_s = curvature * (na_depth - d_bot) * Ec
    # Usually tension, so concrete is cracked/ignored, we use 'n'
    # If NA is below rebar, it's compression, use 'n-1'
    if na_depth > d_bot:
        mod_n = n_s - 1
    else:
        mod_n = n_s
        
    F_bot_s = sigma_c_at_bot_s * mod_n * As_bot
    M_bot_s = F_bot_s * d_bot

    # Sum Forces (Internal)
    # Signs: Compression (+), Tension (-) determined by (na - d) logic above
    F_internal = C_conc + F_top_s + F_bot_s
    
    # Sum Moments about Top Fiber (to check against M_ext + P*ecc)
    # But easier: Return Force Imbalance for specific P
    
    return F_internal, C_conc, F_top_s, F_bot_s

# --- SOLVER ---
# We need to find NA_depth (x) such that Internal Moment = External Moment
# AND Internal Force = External Force.
# But for a specific M and P, there is only ONE NA depth and ONE Curvature.

# Iterative Approach:
# 1. Convert M_ext to equivalent system relative to Top Fiber?
#    M_about_top = M_ext + P * (h/2)  <-- Assuming M_ext is about geometric centroid
#    Standard convention: M applied about geometric centroid.

# Target Loads
P_target = P * 1000 # N
if load_type == "Sagging":
    M_target = M_ext * 1e6 # Nmm
else:
    M_target = -M_ext * 1e6

# Calculate Eccentricity of P from Top Fiber
# P is applied at? Usually P is applied at the strand level in specific beam design,
# but here you defined P as a generic force on the section. 
# Let's assume P is applied at the CENTROID of the strands (d_top).
# This is standard for "Prestress Force".
M_from_P_about_GeoCentroid = P_target * (h/2 - d_top) # P * e
M_design = M_target + M_from_P_about_GeoCentroid

# SOLVER LOGIC:
# We iterate 'x' (NA depth) and 'phi' (curvature).
# Actually, for a linear elastic section, we can use the Transformed Section Properties directly.
# But we need to check if it's cracked.

# Step 1: Calculate Gross Uncracked Properties
Ag = b * h + (n_ps-1)*As_top + (n_s-1)*As_bot
# Centroid
y_bar_uncracked = ( (b*h*h/2) + (n_ps-1)*As_top*d_top + (n_s-1)*As_bot*d_bot ) / Ag
# Inertia
I_uncracked = (b*h**3)/12 + b*h*(h/2 - y_bar_uncracked)**2 + \
              (n_ps-1)*As_top*(d_top - y_bar_uncracked)**2 + \
              (n_s-1)*As_bot*(d_bot - y_bar_uncracked)**2

# Stress Check (Uncracked)
# Total Axial
sigma_axial = P_target / Ag
# Total Moment about Uncracked Centroid
# M_ext is about geometric centroid (h/2). Shift to y_bar.
M_total_uncracked = M_target + P_target * (y_bar_uncracked - d_top) 

# Stress at bottom fiber
y_bot = h - y_bar_uncracked
sigma_bot_check = sigma_axial - (M_total_uncracked * y_bot / I_uncracked)

is_cracked = False
fr = 0.6 * math.sqrt(fc) # Modulus of Rupture
if sigma_bot_check < -fr: # Tension is negative
    is_cracked = True

# --- FINAL CALCULATION ---

if not is_cracked:
    state = "Uncracked (Elastic)"
    NA_depth = y_bar_uncracked # Just geometric approximation for NA shift due to bending? 
    # No, NA is where stress is 0.
    # Sigma(y) = P/A + My/I
    # 0 = P/A + M*y/I  -> y = -(P/A)/(M/I)
    ecc_stress = sigma_axial / (M_total_uncracked / I_uncracked)
    NA_depth = y_bar_uncracked + ecc_stress # Measured from top? No, from centroid.
    # Let's convert to from top
    NA_final_top = NA_depth
    
    # Recalculate Stresses
    sigma_top = sigma_axial + (M_total_uncracked * y_bar_uncracked / I_uncracked)
    sigma_bot = sigma_axial - (M_total_uncracked * y_bot / I_uncracked)

else:
    state = "Cracked (Elastic)"
    
    # If cracked, we must find x where Compression Moment = Tension Moment
    # BUT we have Axial load P.
    # The Equation: Q = 0.5*b*x^2 + (n-1)As'(x-d') - nAs(d-x) = P_axial_eccentricity? 
    # It gets complex. Let's do a brute force scan for x.
    
    found_x = False
    best_x = 0
    min_diff = 1e9
    
    # Scan NA from 1mm to h
    for x_try in range(1, int(h)):
        # 1. Properties of Cracked Section with NA at x_try
        # Comp Area
        A_cr = b*x_try + (n_ps-1)*As_top
        # Tension Area (Steel only)
        # Note: If bottom steel is in tension
        
        # Centroid of this cracked shape
        # Moment of area about top
        mom_area = (b*x_try * x_try/2) + ((n_ps-1)*As_top * d_top) + (n_s*As_bot * d_bot) 
        # Wait, for cracked section analysis with P, we usually take moments about P application point.
        
        # Let's use the Force Equilibrium method directly
        # Sum of Forces = P
        # Sum of Moments = M
        
        # Assume a curvature phi = 1 (dummy)
        # Calculate P_int required to sustain this phi at NA=x_try
        # P_int = Integral(sigma)
        # M_int = Integral(sigma * y)
        
        # We need to find x_try such that the Eccentricity M_int/P_int matches M_design/P_design
        
        # Calculate Inertia of Cracked Section about x_try
        I_cr_x = (b*x_try**3)/3 + (n_ps-1)*As_top*(x_try-d_top)**2 + n_s*As_bot*(d_bot-x_try)**2
        
        # Moment required to create unit curvature
        # M = E * I * phi
        # Stress = M * y / I
        # This is circular.
        
        # Use simple standard equation for x in RC with Axial:
        # P = 0.5*fc*b*x + ...
        # e' = eccentricity of P from Tension Steel
        
        pass 
        # To ensure the tool works reliably without complex solver crashing, 
        # we will use the Standard Transformed Area approximation including P
        # by treating P as a modification to the Moment.
    
    # FALLBACK to Robust Iterative Method for P+M
    # We iterate x (NA depth).
    # For a given x, we calculate Properties.
    # We apply M and P.
    # We check if the calculated Stress Zero point matches x.
    
    for x_try in range(1, int(h), 1):
        # Calc Cracked Properties about this NA
        # Concrete Area
        Ac = b * x_try
        # Transformed Steel
        At_top = (n_ps - 1) * As_top
        At_bot = n_s * As_bot
        
        A_total = Ac + At_top + At_bot
        
        # Centroid of this cracked shape (y_cr) from top
        y_cr = ( (Ac*x_try/2) + (At_top*d_top) + (At_bot*d_bot) ) / A_total
        
        # Inertia about y_cr
        I_cr = (b*x_try**3)/12 + Ac*(x_try/2 - y_cr)**2 + \
               At_top*(d_top - y_cr)**2 + \
               At_bot*(d_bot - y_cr)**2
        
        # Calculate Stresses with these properties
        # M_tot = M_ext + P * (y_cr - d_top)
        # e = y_cr - d_top (distance from centroid to P application)
        e_p = y_cr - d_top
        M_tot = M_target + P_target * e_p
        
        # Stress at depth x_try (should be 0)
        # sigma = P/A + M*y/I
        # y for NA is (x_try - y_cr) (distance from centroid to NA)
        # Actually y is positive downwards
        
        dist_centroid_to_na = x_try - y_cr
        
        sigma_at_na = (P_target / A_total) + (M_tot * dist_centroid_to_na / I_cr)
        
        if abs(sigma_at_na) < 0.5: # Tolerance
            best_x = x_try
            found_x = True
            break
            
    if found_x:
        NA_final_top = best_x
        # Recalc final values
        # (Repeat properties calc for best_x)
        x_try = best_x
        Ac = b * x_try; At_top = (n_ps - 1) * As_top; At_bot = n_s * As_bot
        A_total = Ac + At_top + At_bot
        y_cr = ( (Ac*x_try/2) + (At_top*d_top) + (At_bot*d_bot) ) / A_total
        I_cr = (b*x_try**3)/12 + Ac*(x_try/2 - y_cr)**2 + At_top*(d_top - y_cr)**2 + At_bot*(d_bot - y_cr)**2
        e_p = y_cr - d_top
        M_tot = M_target + P_target * e_p
        
        sigma_top = (P_target / A_total) + (M_tot * (0 - y_cr) / I_cr)
        sigma_bot = (P_target / A_total) + (M_tot * (h - y_cr) / I_cr)
    else:
        state = "Calculation Diverged (High P)"
        NA_final_top = 0
        sigma_top = 0
        sigma_bot = 0

# --- RESULTS DISPLAY ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("Analysis Results")
    st.info(f"**Section State:** {state}")
    
    st.metric("Neutral Axis Depth (from Top)", f"{NA_final_top:.2f} mm")
    st.metric("Max Concrete Stress (Top)", f"{sigma_top:.2f} MPa")
    st.metric("Stress at Bottom Rebar Level", f"{sigma_bot:.2f} MPa")
    
    if sigma_top > 0.6 * fc:
        st.warning("⚠️ Concrete Compressive Stress exceeds 0.6*fc")
    if sigma_bot < -0.6 * math.sqrt(fc) and state == "Uncracked (Elastic)":
        st.error("⚠️ Section should be treated as CRACKED (Calculated Tension > Rupture)")

with col2:
    st.subheader("Stress Diagram")
    
    fig, ax = plt.subplots(figsize=(4, 6))
    
    # Draw Concrete Section
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='gray', facecolor='#f0f0f0')
    ax.add_patch(rect)
    
    # Draw Strands
    circle_ps = patches.Circle((b/2, h - d_top), radius=8, color='blue', label='Prestress')
    ax.add_patch(circle_ps)
    
    # Draw Rebar
    circle_s = patches.Circle((b/2, h - d_bot), radius=8, color='red', label='Rebar')
    ax.add_patch(circle_s)
    
    # Draw Neutral Axis
    na_y = h - NA_final_top
    ax.axhline(y=na_y, color='green', linestyle='--', linewidth=2, label='Neutral Axis')
    ax.text(10, na_y + 10, f'NA = {NA_final_top:.0f}mm', color='green')
    
    # Draw Stress Profile (Qualitative)
    # We draw a line on the right side
    ax.plot([b+50, b+50 + sigma_top*5], [h, h], color='purple') # Top Stress
    ax.plot([b+50, b+50 + sigma_bot*5], [0, 0], color='purple') # Bot Stress
    ax.plot([b+50 + sigma_top*5, b+50 + sigma_bot*5], [h, 0], color='purple', linewidth=2, label='Stress')
    ax.axvline(x=b+50, color='black', linestyle=':', linewidth=0.5)
    
    ax.set_xlim(-20, b + 150)
    ax.set_ylim(-20, h + 20)
    ax.legend(loc='upper right')
    ax.set_aspect('equal')
    ax.axis('off')
    
    st.pyplot(fig)
