# How to Compile and Use the UDF in ANSYS Fluent 🛠️

Follow the steps below to compile and apply the soot formation UDF in ANSYS Fluent.

---

## 1. Prepare the Simulation Case

Before compiling the UDF, ensure your Fluent case is properly configured:

- Enable:
  - Energy equation  
  - Species transport  

- Your reaction mechanism must include:
  - C₂H₂ (soot precursor)  
  - H₂  
  - H

---

## 2. Enable User-Defined Scalars (UDS)

The soot model uses two additional transport equations.

- Navigate to:

  Models → User-Defined → Scalars


- Enable **2 UDS equations** and assign:
  - **UDS-0 → Soot number density (N)**  
  - **UDS-1 → Soot carbon concentration (C)**  

---

## 3. Compile (Interpret) the UDF

This implementation uses the **interpreted UDF approach**.

1. Navigate to:
   ```
   User-Defined → Functions → Interpreted
   ```

2. Click **Browse...** and select:
   ```
   AR_version_monodisperse_agglomerate_udf_uds.c
   ```

3. Click **Interpret**

✅ If there are no errors in the console, the UDF has been successfully loaded.

---

## 4. Assign Source Terms (Hook the UDF)

Attach the source terms to the appropriate scalar equations.

- Navigate to:
  ```
  Physics → Cell Zone Conditions → (select your fluid zone) → Edit
  ```

- Go to the **Source Terms** tab

### For UDS-0 (Soot Number Density)
Assign:
- `N_nucl` → nucleation  
- `N_coag` → coagulation  

### For UDS-1 (Soot Carbon Concentration)
Assign:
- `C_nucl` → nucleation  
- `C_SG` → surface growth  

Click **Apply** to confirm.

---

## 5. Set Scalar Diffusivity

- Navigate to:
  ```
  Physics → Materials → Create/Edit
  ```

- Select:
  - **Material Type → Mixture**

- In the **Properties** section:
  - Locate **UDS Diffusivity**
  - Select **User-Defined**
  - Click **Edit...**
  - Assign:
    ```
    scalar_diff
    ```

---

## 6. Boundary Conditions, Initialization, and Solution

### Check Boundary Conditions
- Ensure UDS variables are properly defined at:
  - Inlets  
  - Walls (if required)  

---

### Initialize the Solution

- Navigate to:
  ```
  Solution → Initialization
  ```
- Click **Initialize**

---

### Run the Simulation

---
