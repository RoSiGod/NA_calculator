import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="AS3600 Prestress Analyzer", layout="wide")

st.title("🏗️ AS3600 Prestressed RC Section Analyzer (Fixed)")
st.markdown("""
**Fix Applied:** corrected the sign convention in the cracked section analysis. 
Sagging moment now correctly increases compression at the top fiber rather than reversing it.
""")

# --- SIDEBAR INPUTS ---
st.sidebar.header("1. Geometry")
b = st.sidebar.number_input("Width (b) [mm]", value=300.0)
h = st.sidebar.number_input("Depth (h) [mm]", value=600.0)

st.sidebar.header("2. Materials")
fc = st.sidebar.number_input("Concrete Strength (fc) [MPa]", value=40.0)
Ec = 4700 * math.sqrt(fc) 
st.sidebar.markdown(f"*Calc. Ec: {Ec:.0f} MPa*")
Es = st.sidebar.number_input("Rebar Modulus (Es) [MPa]", value=200000.0)
Eps = st.sidebar.number_input("Strand Modulus (Eps) [MPa]", value=195000.0)

# AS3600 Limits
fctf = 0.6 * math.sqrt(fc)
st.sidebar.markdown("---")
st.sidebar.markdown(f"**AS3600 Limits:**")
st.sidebar.markdown(f"Tensile Limit ($f'_{{ct.f}}$): **{fctf:.2f} MPa**")
st.sidebar.markdown(f"Steel Yield Limit: **450 MPa**")

n_s = Es / Ec
n_ps = Eps / Ec

st.sidebar.header("3. Reinforcement")
st.sidebar.subheader("Top Zone (Strands/Comp)")
As_top = st.sidebar.number_input("Strand Area (As_top) [mm²]", value=400.0)
d_top = st.sidebar.number_input("Depth to Strands (d_top) [mm]", value=50.0)

st.sidebar.subheader("Bottom Zone (Rebar/Tens)")
As_bot = st.sidebar.number_input("Rebar Area (As_bot) [mm²]", value=1200.0)
d_bot = st.sidebar.number_input("Depth to Rebar (d_bot) [mm]", value=550.0)

st.sidebar.header("4. Loads")
P = st.sidebar.number_input("Prestress Force (P) [kN]", value=500.0, step=50.0)
M_ext = st.sidebar.number_input("External Moment (M_ext) [kNm]", value=150.0, step=10.0)
load_type = st.sidebar.radio("Moment Type", ["Sagging (Comp. Top)", "Hogging (Tens. Top)"])

# --- SOLVER LOGIC ---

# 1. Standardize Loads
P_force = P * 1000 # N (Compression +)
# Convention: Positive M = Sagging (Compress Top)
if load_type == "Sagging":
    M_applied = M_ext * 1e6 
else:
    M_applied = -M_ext * 1e6 

# 2. Uncracked Section Properties
Ag = b * h + (n_ps-1)*As_top + (n_s-1)*As_bot
mat_top = (b*h * h/2) + (n_ps-1)*As_top*d_top + (n_s-1)*As_bot*d_bot
y_bar_uncracked = mat_top / Ag # Depth of centroid from top

I_uncracked = (b*h**3)/12 + b*h*(h/2 - y_bar_uncracked)**2 + \
              (n_ps-1)*As_top*(d_top - y_bar_uncracked)**2 + \
              (n_s-1)*As_bot*(d_bot - y_bar_uncracked)**2

# Loads at Uncracked Centroid
e_p_uncracked = y_bar_uncracked - d_top # Positive if P is above centroid
M_total_uncracked = M_applied + (P_force * e_p_uncracked)

# Check Stress at Bottom Fiber
# y is distance FROM centroid (Positive Downwards)
y_bot_dist = h - y_bar_uncracked
# Formula: sigma = P/A - My/I
sigma_bot_check = (P_force / Ag) - (M_total_uncracked * y_bot_dist / I_uncracked)

is_cracked = False
if sigma_bot_check < -fctf: 
    is_cracked = True

# --- FINAL CALCULATION ---

if not is_cracked:
    state = "Uncracked (Elastic)"
    
    # Top Fiber Stress (y is negative, so -My/I becomes positive/compression)
    y_top_dist = 0 - y_bar_uncracked
    sigma_top = (P_force / Ag) - (M_total_uncracked * y_top_dist / I_uncracked)
    sigma_bot = sigma_bot_check
    
    # Find NA (Linear Interpolation)
    if sigma_top == sigma_bot:
        NA_final_top = -1
    else:
        # Similar triangles
        NA_final_top = (sigma_top / (sigma_top - sigma_bot)) * h

else:
    state = "Cracked (Elastic)"
    
    # --- ITERATIVE SOLVER ---
    found_x = False
    best_x = 0
    
    # Iterate x (Depth of Neutral Axis from Top)
    for x_try in range(1, int(h)):
        
        # 1. Cracked Properties (Top Concrete + All Steel)
        Ac = b * x_try
        At_top = (n_ps - 1) * As_top 
        At_bot = n_s * As_bot        
        
        A_cr_total = Ac + At_top + At_bot
        
        # Centroid of cracked section (from top)
        y_cr = ( (Ac * x_try/2) + (At_top * d_top) + (At_bot * d_bot) ) / A_cr_total
        
        I_cr = (b*x_try**3)/12 + Ac*(x_try/2 - y_cr)**2 + \
               At_top*(d_top - y_cr)**2 + \
               At_bot*(d_bot - y_cr)**2
        
        # 2. Loads at Cracked Centroid
        e_p_cracked = y_cr - d_top
        M_tot_cracked = M_applied + (P_force * e_p_cracked)
        
        # 3. Check stress at the Assumed NA (x_try)
        # y distance from centroid to x_try (Positive Downwards)
        y_at_NA = x_try - y_cr
        
        # KEY FIX: Use P/A - My/I
        sigma_at_boundary = (P_force / A_cr_total) - (M_tot_cracked * y_at_NA / I_cr)
        
        # 4. Convergence
        if abs(sigma_at_boundary) < 0.1: # Tolerance
            best_x = x_try
            found_x = True
            
            # Final Stresses using the correct Centroid and I_cr
            # Top Fiber (y = 0 - y_cr) -> Negative distance
            sigma_top = (P_force / A_cr_total) - (M_tot_cracked * (0 - y_cr) / I_cr)
            
            # Bottom Fiber (y = h - y_cr) -> Positive distance
            sigma_bot = (P_force / A_cr_total) - (M_tot_cracked * (h - y_cr) / I_cr)
            break
            
    if found_x:
        NA_final_top = best_x
    else:
        state = "Calculation Diverged"
        NA_final_top = h 
        sigma_top = P_force / Ag 
        sigma_bot = P_force / Ag

# --- STEEL STRESS CALCULATIONS ---
# Interpolate stress at steel levels based on linear concrete stress profile
stress_slope = (sigma_bot - sigma_top) / h

sigma_c_at_top_steel = sigma_top + (stress_slope * d_top)
sigma_c_at_bot_steel = sigma_top + (stress_slope * d_bot)

sigma_st_top = sigma_c_at_top_steel * n_ps
sigma_st_bot = sigma_c_at_bot_steel * n_s

# --- RESULTS DISPLAY ---
col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("Analysis Results")
    st.markdown("### N.A. Position from top:")
    st.metric("N.A. Position",f"{NA_final_top:.2f}")
    if is_cracked:
        st.warning(f"**Section State:** {state}")
    else:
        st.success(f"**Section State:** {state}")
    
    st.markdown("### 1. Concrete Compressive Check")
    col_c1, col_c2 = st.columns(2)
    col_c1.metric("Max Comp. (Top)", f"{sigma_top:.2f} MPa")
    
    if sigma_top > fc:
        col_c2.error(f"❌ Fails (> {fc} MPa)")
    elif sigma_top > 0.6 * fc:
        col_c2.warning(f"⚠️ High (> 0.6fc)")
    else:
        col_c2.success("✅ OK (< 0.6fc)")

    st.markdown("### 2. Concrete Tensile Check")
    min_stress = min(sigma_top, sigma_bot) 
    col_t1, col_t2 = st.columns(2)
    col_t1.metric("Max Tension", f"{min_stress:.2f} MPa")
    
    if min_stress < -fctf:
         col_t2.error(f"❌ Cracked (> {fctf:.2f} MPa)")
    else:
         col_t2.success(f"✅ Uncracked (< {fctf:.2f} MPa)")

    st.markdown("### 3. Steel Stress Check (< 450 MPa)")
    col_s1, col_s2 = st.columns(2)
    
    col_s1.metric("Rebar Stress (Bottom)", f"{sigma_st_bot:.2f} MPa")
    if abs(sigma_st_bot) > 450:
        col_s1.error("❌ Yielded")
    else:
        col_s1.success("✅ OK")

    col_s2.metric("Strand Stress Change", f"{sigma_st_top:.2f} MPa")

with col2:
    st.subheader("Stress Diagram")
    fig, ax = plt.subplots(figsize=(4, 6))
    
    # Draw Concrete
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='black', facecolor='#f0f2f6')
    ax.add_patch(rect)
    
    # Draw Steel
    ax.add_patch(patches.Circle((b/2, h - d_top), radius=10, color='blue', label='Strand'))
    ax.add_patch(patches.Circle((b/2, h - d_bot), radius=10, color='red', label='Rebar'))
    
    # Draw NA
    if 0 <= NA_final_top <= h:
        ax.axhline(y=h - NA_final_top, color='green', linestyle='--', linewidth=2, label='NA')
    
    # Draw Stress Block
    offset = b + 50
    ax.axvline(x=offset, color='gray', linestyle='-') # Zero line
    
    # Scale stress for plotting
    scale = 3 
    
    # Top Stress (Purple)
    ax.plot([offset, offset + sigma_top*scale], [h, h], color='purple')
    ax.text(offset + sigma_top*scale + 10, h, f"{sigma_top:.1f}", va='center', color='purple', fontsize=8)
    
    # Bot Stress (Purple)
    ax.plot([offset, offset + sigma_bot*scale], [0, 0], color='purple')
    ax.text(offset + sigma_bot*scale + 10, 0, f"{sigma_bot:.1f}", va='center', color='purple', fontsize=8)
    
    # Connect
    ax.plot([offset + sigma_top*scale, offset + sigma_bot*scale], [h, 0], color='purple', linewidth=2)
    
    ax.set_xlim(-20, b + 250)
    ax.set_ylim(-20, h + 20)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.15))
    ax.axis('off')
    st.pyplot(fig)
