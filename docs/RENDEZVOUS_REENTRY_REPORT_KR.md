# Rendezvous and Re-entry Mission Simulator 보고서

## 요약

본 코드는 Earth orbit의 unmanned chaser가 target spacecraft와 rendezvous / proximity operation을 수행한 뒤, 120 km entry interface에서 re-entry vehicle을 분리하여 atmospheric re-entry를 해석하는 mission-level simulator이다. MATLAB은 전체 mission sequence, 3-DOF/6-DOF propagation, proximity operation, de-orbit, re-entry diagnostics를 담당하고, Python은 Phase 1 `CUSTOM_IMPULSE`에 들어가는 `phase angle`, `delta-V`, `gamma`를 shooting / optimization 방식으로 탐색한 뒤 MATLAB-readable JSON으로 넘기는 역할을 한다.

현재 active code의 핵심 선택지는 다음과 같다.

- Phase 1 phasing mode: `CUSTOM_IMPULSE`, `HOHMANN`, `MULTI_HOHMANN`
- Burn model: `IMPULSIVE`, `FINITE_BURN`
- Orbital atmospheric drag: `OFF`, `ISA76`
- Phase 3 re-entry setup: `HOHMANN`, `R_BAR_200_FPA`
- Phase 4 re-entry vehicle shape: `COMPROMISE`, `HEATLOAD_MIN`, `PAYLOAD_MAX`, `TPS_MIN`

방치되어 있던 `LAMBERT` mode는 active code에서 제거했다. 또한 과거의 `continuous` burn alias는 실제 low-thrust spiral 모델과 혼동될 수 있으므로 더 이상 허용하지 않는다. 짧은 시간 동안 impulsive delta-V를 실제 thrust로 분산하는 옵션은 명확하게 `FINITE_BURN`으로만 사용한다.

---

# Part I. 이론적 배경

## 1. Mission phase 개요

전체 mission은 네 개의 phase로 나뉜다.

1. **Phase 1: Phasing & Homing**
   Chaser가 insertion orbit에서 target 근처의 waiting / capture geometry로 접근한다. Python shooting으로 얻은 `phase angle`, `delta-V`, `gamma`를 MATLAB의 `CUSTOM_IMPULSE` mode에 넣는 것이 현재 권장 workflow이다.

2. **Phase 2: Proximity Operation**
   Target-centered `LVLH frame`에서 S2, S3, S4 waypoint를 따라 접근한다. 기본 구조는 S2 hold trim, cycloidal drift, R-bar hop, braking impulse로 구성된다.

3. **Phase 3: De-orbit / Entry Interface Setup**
   Vehicle을 120 km atmospheric entry interface로 보낸다. `HOHMANN`은 direct FPA-targeted descent를 수행하고, `R_BAR_200_FPA`는 200 km parking orbit과 R-bar alignment를 거친 뒤 120 km interface로 진입한다.

4. **Phase 4: Atmospheric Entry**
   120 km부터 re-entry vehicle을 별도 객체로 전파한다. 동시에 chaser는 orbiting relay로 계속 propagation되며, re-entry vehicle과 chaser 사이의 geometric line of sight가 유지되는지 확인한다.

## 2. Coordinate frame

### 2.1 ECI-like inertial frame

MATLAB propagation의 기본 position / velocity state는 Earth-centered inertial에 가까운 frame으로 취급한다. Central gravity, J2 perturbation, orbital propagation, Phase 1/3 transfer는 이 frame에서 계산된다.

State vector는 상황에 따라 두 가지 형태로 사용된다.

- 14-state chaser:
  `[r(3); v(3); q(4); w(3); mass]`
- 6-state target:
  `[r(3); v(3)]`

여기서 `r`은 ECI position, `v`는 ECI velocity, `q`는 quaternion, `w`는 angular velocity이다.

### 2.2 LVLH frame

Proximity operation은 target-centered `LVLH (Local Vertical Local Horizontal)` frame에서 해석한다.

- `R-bar`: target radial direction
- `V-bar`: along-track / velocity direction
- `H-bar`: orbital angular momentum direction

본 코드의 Phase 2는 target state로부터 LVLH basis를 만들고, chaser-target relative position / velocity를 해당 frame으로 변환한다. Relative velocity는 단순히 `v_chaser - v_target`만 쓰지 않고, LVLH frame rotation effect를 포함한다.

### 2.3 Co-rotating atmosphere

Atmospheric drag와 re-entry aerodynamics에서는 atmosphere-relative velocity가 중요하다. 코드에서는 전체 state를 ECI로 유지하되, aerodynamic velocity만 다음처럼 계산한다.

```text
v_rel = v_ECI - omega_Earth x r
```

즉, full ECEF state transition을 새로 도입하지 않고, co-rotating atmosphere에 대한 상대속도를 사용한다. 현재 fidelity에서는 가장 깔끔한 compromise이다.

## 3. Orbital dynamics

### 3.1 Central gravity

기본 translational acceleration은 two-body gravity이다.

```text
a_g = -mu / r^3 * r
```

여기서 `mu`는 Earth gravitational parameter, `r`은 geocentric distance이다.

### 3.2 J2 perturbation

Earth oblateness는 가장 큰 zonal harmonic perturbation인 `J2`로 모델링한다. 코드의 `Env_EOM.m`은 다음 형태의 Cartesian J2 acceleration을 더한다.

```text
a_J2 = 1.5 * J2 * mu/r^2 * (Re/r)^2 *
       [ x/r * (5z^2/r^2 - 1);
         y/r * (5z^2/r^2 - 1);
         z/r * (5z^2/r^2 - 3) ]
```

J2는 LEO에서 RAAN, argument of perigee, phase evolution에 누적 오차를 만든다. 그래서 Phase 1/3의 wait-time search는 단순 Keplerian closed-form만 쓰지 않고, chaser와 target을 같이 numerical propagation하면서 arrival geometry를 찾는다.

## 4. Maneuver model

### 4.1 Impulsive burn

`IMPULSIVE` mode는 delta-V가 순간적으로 적용된다고 가정한다.

```text
v_plus = v_minus + delta_v_vector
```

Mass depletion은 ideal rocket equation으로 계산한다.

```text
m_f = m_0 * exp(-DeltaV / (Isp * g0))
```

이 모델은 optimization과 mission-level budget 계산에 빠르고 안정적이다. 하지만 실제 engine firing time, finite burn gravity loss, attitude coupling은 생략된다.

### 4.2 Finite burn

`FINITE_BURN` mode는 같은 delta-V vector를 짧은 burn duration 동안 fixed direction thrust로 나누어 적용한다. 이것은 low-thrust spiral이 아니다. 사용자가 의도한 것처럼 impulsive maneuver를 실제 high-thrust firing duration으로 펼친 모델이다.

Burn duration은 mass flow와 rocket equation으로 근사한다.

```text
ve = Isp * g0
m_f = m_0 * exp(-DeltaV / ve)
burn_time = (m_0 - m_f) * ve / thrust
```

현재 finite burn은 maneuver direction을 burn 시작 시점에 고정한다. 아주 긴 burn이나 attitude steering이 중요한 경우에는 향후 guidance-coupled finite burn model이 필요하다.

### 4.3 Multi-Hohmann과 thermal constraint

긴 burn은 chamber overheating이나 actuator thermal limit을 유발할 수 있다. `MULTI_HOHMANN` mode는 큰 altitude change를 여러 Hohmann leg로 분할한다.

코드는 다음 조건을 만족하는 leg count를 자동으로 찾을 수 있다.

- maximum single-burn duration
- maximum single-burn delta-V

따라서 같은 total transfer라도 single burn peak demand를 낮출 수 있다.

## 5. Hohmann transfer와 phase angle

Two-body Hohmann transfer의 half-period transfer time은 다음과 같다.

```text
a_t = (r1 + r2) / 2
TOF = pi * sqrt(a_t^3 / mu)
```

Target이 circular orbit에 있다고 보면, transfer 중 target이 이동하는 angle은 `n2 * TOF`이고, chaser가 opposite side에 도달해야 하므로 required phase angle은 다음과 같이 쓸 수 있다.

```text
phase = pi - n2 * TOF
      = pi * (1 - (a_t / r2)^(3/2))
```

`Mission_Config.m`의 `sys.phase`는 이제 SI 단위 일관성을 유지하여 이 식으로 계산된다. 300 km에서 495 km로 올라가는 기본 설정에서는 약 3.8165 deg이다.

다만 실제 code는 J2와 target co-propagation을 포함하므로, closed-form phase angle은 initial guess / filter에 가깝고 최종 alignment는 numerical search가 담당한다.

## 6. Phase 2 relative motion

Phase 2는 `LVLH frame`에서 waypoint impulse 접근을 한다. 기본 철학은 상태를 순간적으로 덮어쓰지 않고, impulse와 free propagation만으로 상태를 이동시키는 것이다.

### 6.1 CW / HCW targeting

Target orbit이 circular이고 perturbation이 작다고 보면, local relative motion은 `Clohessy-Wiltshire` 또는 `Hill-Clohessy-Wiltshire (HCW)` equation으로 근사할 수 있다.

Planar form의 대표적인 구조는 다음과 같다.

```text
x_ddot - 2n y_dot - 3n^2 x = 0
y_ddot + 2n x_dot = 0
z_ddot + n^2 z = 0
```

코드는 CW state transition relation을 이용해 waypoint까지 필요한 departure delta-V를 계산하고, 실제 propagation은 nonlinear `Env_EOM.m`으로 수행한다. 즉, CW는 guidance guess이고, 실제 motion은 J2 포함 nonlinear dynamics로 검증된다.

### 6.2 Cycloidal drift

S2에서 V-bar 방향으로 작은 impulse를 주면 natural relative motion에 의해 R-bar excursion이 생긴다. 코드에서는 다음 관계를 사용한다.

```text
|Delta_R|max = 4 * v0 / n
```

이후 V-bar coordinate가 0을 처음 crossing하는 지점을 S3로 잡고, 거기서 residual relative velocity를 braking impulse로 제거한다.

### 6.3 R-bar hop

S3부터 S4까지는 여러 R-bar waypoint로 나누어 이동한다. 각 hop은 다음 구조를 가진다.

1. 현재 LVLH relative state 계산
2. CW targeting으로 departure impulse 계산
3. nonlinear propagation
4. arrival residual velocity braking

이 구조는 하나의 큰 approach보다 더 해석하기 쉽고, waypoint별 arrival error와 braking delta-V를 확인할 수 있다.

## 7. Phase 3 entry interface design

Phase 3의 목적은 re-entry vehicle을 120 km atmospheric entry interface에 보내는 것이다. 여기서 중요한 값은 altitude뿐 아니라 `flight path angle (FPA)`이다.

```text
FPA = atan2(v_radial, v_horizontal)
```

현재 target FPA magnitude는 `Mission_Config.m`의 `sys.reentry_flight_path_angle`로 설정되며 기본값은 4 deg이다. Descending entry이므로 실제 출력 FPA는 약 -4 deg가 된다.

### 7.1 HOHMANN mode

`HOHMANN` Phase 3 mode는 현재 상태에서 120 km / target FPA interface로 향하는 de-orbit injection을 수행하고, 120 km crossing에서 propagation을 멈춘다. 예전처럼 120 km 아래로 내려간 뒤 nonphysical circularization을 수행하지 않는다.

### 7.2 R_BAR_200_FPA mode

`R_BAR_200_FPA` mode는 세 단계로 구성된다.

1. 200 km parking orbit으로 lowering
2. target 기준 R-bar 아래쪽 alignment까지 coast
3. 120 km / target FPA interface injection

이 mode는 re-entry 직전 geometry를 target-relative 관점에서 더 명확히 만들기 위한 option이다.

## 8. Atmospheric drag

Drag acceleration은 다음 식을 사용한다.

```text
a_D = -0.5 * rho * Cd * A / m * |v_rel| * v_rel
```

여기서 `rho`는 atmosphere density, `Cd`는 drag coefficient, `A`는 reference area, `m`은 mass, `v_rel`은 co-rotating atmosphere-relative velocity이다.

### 8.1 ISA76 atmosphere

MATLAB에서 `atmosisa`가 사용 가능하고 altitude가 lower atmosphere 범위에 있으면 이를 우선 사용한다. 그 외에는 `Standard_Atmosphere_Density.m`의 ISA76-style fallback이 사용된다.

120 km 이상에서는 density가 매우 작지만, re-entry가 진행되면서 dynamic pressure와 heating이 급격히 증가한다.

## 9. Re-entry aerodynamics

Re-entry vehicle shape는 사용자가 제공한 PPT의 geometry와 L/D graph를 기반으로 한다. 현재 code에 들어간 shape는 다음 네 가지이다.

- `COMPROMISE`
- `HEATLOAD_MIN`
- `PAYLOAD_MAX`
- `TPS_MIN`

각 shape는 length, max width, max height, base diameter, nose radius, aspect ratio, reference area, approximate L/D curve를 가진다.

현재 limitation은 `Cd`가 shape별 high-fidelity aerodynamic table이 아니라 constant assumption이라는 점이다. PPT에는 geometry와 L/D trend는 있지만 Mach/AoA별 `Cd`, `Cl`, `Cm` table은 없으므로, 지금 모델은 trajectory-level sensitivity용으로 보는 것이 맞다.

Lift는 drag magnitude와 L/D lookup으로 계산된다.

```text
L = D * (L/D)
```

Lift direction은 local vertical plane을 기준으로 만들고, `bank_angle_deg`로 회전시킬 수 있다. 기본 bank angle은 0 deg이다.

## 10. Heating, dynamic pressure, g-load

### 10.1 Dynamic pressure

Dynamic pressure는 다음과 같다.

```text
q = 0.5 * rho * V_rel^2
```

Max dynamic pressure는 구조하중과 control authority를 판단하는 첫 번째 지표이다.

### 10.2 Sutton-Graves heat flux

Stagnation-point convective heat flux는 Sutton-Graves 형태의 engineering correlation으로 계산한다.

```text
q_dot = k * sqrt(rho / R_n) * V_rel^3
```

여기서 `R_n`은 nose radius, `k`는 Earth entry용 coefficient이다. 코드 기본값은 `1.83e-4`이며 SI unit 기준으로 사용한다.

Heat load는 heat flux를 time integration한 값이다.

```text
Q = integral(q_dot dt)
```

이는 TPS sizing의 첫 번째 order estimate로 유용하지만, radiation heating, ablation, wall catalysis, real-gas nonequilibrium, shock-layer chemistry는 포함하지 않는다.

### 10.3 Aero g-load

Aero g-load는 aerodynamic acceleration magnitude를 `g0`로 나눈 값이다.

```text
g_load = |a_drag + a_lift| / g0
```

Max g-load는 payload survivability와 structure sizing에서 중요한 diagnostic이다.

## 11. Line-of-sight communication check

Phase 4에서 chaser는 re-entry vehicle과 분리된 뒤에도 orbiting relay로 계속 propagation된다. Re-entry vehicle은 atmosphere 안으로 내려가고, chaser는 orbit에 남는다.

LOS 판단은 두 vehicle을 잇는 line segment와 Earth sphere의 관계로 계산한다.

1. Re-entry vehicle에서 chaser까지 vector를 만든다.
2. 그 segment에서 Earth center에 가장 가까운 점을 찾는다.
3. 그 closest point radius가 `Re + margin`보다 크면 geometric LOS가 clear하다고 본다.

출력 diagnostic은 다음과 같다.

- `los_maintained`: 전체 re-entry 동안 LOS 유지 여부
- `min_los_clearance_m`: Earth limb clearance 최소값
- `min_los_elevation_deg`: re-entry vehicle 기준 chaser elevation 최소값
- `first_los_loss_time_s`: 처음 LOS가 끊긴 시간

중요한 limitation은 plasma blackout physics가 아직 없다는 점이다. 즉, 현재 LOS는 geometry-only check이다. 실제 blackout-free communication을 보려면 electron density, plasma frequency, RF frequency, antenna pointing, link budget까지 추가해야 한다.

---

# Part II. 코드 사용법 및 특징

## 1. 가장 기본 실행

MATLAB에서 repository root를 current folder로 설정한 뒤 실행한다.

```matlab
Main_Mission_Simulator
```

기본값은 `Mission_Config.m`과 `Main_Mission_Simulator.m` 내부 default를 사용한다. 단, Phase 1의 hard-coded `CUSTOM_IMPULSE` 값은 연구 중간값이므로, 정밀한 run은 Python optimizer로 JSON을 만든 뒤 MATLAB에서 불러오는 workflow를 권장한다.

## 2. Python optimizer 실행

Python dependency를 설치한다.

```bash
pip install -r requirements.txt
```

기본 impulsive optimization을 실행하고 MATLAB config를 저장한다.

```bash
python J2PolarHohmannShooting.py --no-plot --matlab-config-out configs/latest_python_solution.json
```

Finite burn까지 포함해서 optimization하려면 다음처럼 실행한다.

```bash
python J2PolarHohmannShooting.py --no-plot --burn-model finite_burn --matlab-config-out configs/latest_python_solution.json
```

Atmospheric drag를 Python optimization에도 반영하고 MATLAB config로 넘기려면 다음처럼 실행한다.

```bash
python J2PolarHohmannShooting.py --no-plot --atmospheric-drag isa76 --matlab-config-out configs/latest_python_solution.json
```

Python script는 output path의 parent directory가 없으면 자동으로 생성한다.

## 3. MATLAB에서 Python JSON 불러오기

MATLAB에서 environment variable을 설정한다.

```matlab
setenv('RENDEZVOUS_CONFIG_JSON','configs/latest_python_solution.json')
Main_Mission_Simulator
```

JSON에 들어있는 값만 override된다. JSON에 없는 값은 `Mission_Config.m`과 script default가 유지된다.

## 4. Burn model 선택

MATLAB에서 burn model을 강제로 바꾸려면 다음을 사용한다.

```matlab
setenv('RENDEZVOUS_BURN_MODEL','IMPULSIVE')
Main_Mission_Simulator
```

또는:

```matlab
setenv('RENDEZVOUS_BURN_MODEL','FINITE_BURN')
Main_Mission_Simulator
```

`CONTINUOUS`는 더 이상 지원하지 않는다. 이 프로젝트에서 finite burn은 low-thrust continuous ascent가 아니라, impulsive delta-V를 수 초~수십 초의 실제 thrust firing으로 펼친 model이다.

## 5. Phase 3 mode 선택

Direct entry interface descent:

```matlab
setenv('RENDEZVOUS_PHASE3_MODE','HOHMANN')
Main_Mission_Simulator
```

200 km parking orbit과 R-bar alignment를 포함하는 mode:

```matlab
setenv('RENDEZVOUS_PHASE3_MODE','R_BAR_200_FPA')
Main_Mission_Simulator
```

`R_BAR_200_FPA`에서 final 200 km -> 120 km injection의 propellant를 budget에 포함하려면:

```matlab
setenv('RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL','on')
Main_Mission_Simulator
```

기본값은 fuel을 제외하되 delta-V는 보고한다.

## 6. Atmospheric drag 선택

Orbital phase에도 ISA76 drag를 켜려면:

```matlab
setenv('RENDEZVOUS_ATMOSPHERIC_DRAG','ISA76')
Main_Mission_Simulator
```

끄려면:

```matlab
setenv('RENDEZVOUS_ATMOSPHERIC_DRAG','OFF')
Main_Mission_Simulator
```

Phase 4 re-entry는 항상 atmosphere model을 사용한다. `RENDEZVOUS_ATMOSPHERIC_DRAG`는 주로 orbital propagation drag option이다.

## 7. Re-entry vehicle shape 선택

기본 shape는 `COMPROMISE`이다. 다른 shape를 사용하려면:

```matlab
setenv('RENDEZVOUS_REENTRY_SHAPE','TPS_MIN')
Main_Mission_Simulator
```

가능한 값:

- `COMPROMISE`
- `HEATLOAD_MIN`
- `PAYLOAD_MAX`
- `TPS_MIN`

## 8. JSON config 예시

아래 JSON은 Python output이 아니어도 직접 작성해서 사용할 수 있다.

```json
{
  "scenario": {
    "initial_chaser_altitude_m": 300000,
    "target_altitude_m": 500000,
    "waiting_altitude_m": 495000,
    "initial_chaser_angle_deg": 0,
    "initial_phase_angle_deg": 90
  },
  "phase1": {
    "mode": "CUSTOM_IMPULSE",
    "burn_model": "FINITE_BURN",
    "phase_angle_deg": 3.8,
    "delta_v_m_s": 58.0,
    "gamma_deg": 0.0
  },
  "maneuver": {
    "finite_burn_thrust_N": 300,
    "finite_burn_isp_s": 200,
    "finite_burn_dt_s": 0.1,
    "max_single_burn_duration_s": 120
  },
  "environment": {
    "atmospheric_drag": {
      "enabled": false,
      "model": "ISA76",
      "chaser_cd": 2.2,
      "chaser_area_m2": 4.0,
      "target_cd": 2.2,
      "target_area_m2": 4.0
    }
  },
  "reentry": {
    "shape": "COMPROMISE",
    "dt_s": 0.5,
    "max_time_s": 2500,
    "terminal_altitude_m": 20000,
    "lift_enabled": true,
    "bank_angle_deg": 0,
    "los_margin_altitude_m": 0
  }
}
```

## 9. Output 해석

MATLAB run이 끝나면 다음이 출력된다.

- Phase별 delta-V
- Phase별 propellant consumption
- Remaining mass
- Phase 1 max burn duration / max single-burn delta-V
- Phase 3 final injection delta-V
- Phase 4 max heat flux
- Phase 4 total heat load
- Phase 4 max dynamic pressure
- Phase 4 max aero g-load
- LOS maintained 여부
- first LOS loss time

주요 figure:

- Proximity Operations: LVLH R-bar / V-bar 접근 궤적
- Mission Comprehensive Dashboard:
  - 3D ECI trajectory
  - altitude history
  - mass history
- Atmospheric Re-entry Diagnostics:
  - altitude vs time
  - atmosphere-relative speed
  - dynamic pressure
  - heat flux
  - aero g-load
  - LOS clearance / elevation

## 10. 주요 파일

### `Main_Mission_Simulator.m`

전체 mission script이다. Phase 1부터 Phase 4까지 순서대로 실행하고, budget table과 plot을 만든다.

### `Mission_Config.m`

Earth constants, orbit altitude, propulsion, maneuver defaults, atmosphere, re-entry vehicle shape, noise placeholder를 정의한다.

### `Phasing_Propagator.m`

Phase 1 phasing / Hohmann / Multi-Hohmann / custom impulse logic을 담당한다. Active code에서 Lambert branch는 제거되었다.

### `Env_EOM.m`

Gravity, J2, optional drag, 6-DOF thrust/mass dynamics를 계산한다.

### `Atmospheric_Drag_Acceleration.m`

Orbital phase drag acceleration을 계산한다.

### `Standard_Atmosphere_Density.m`

ISA76-style density helper이다. MATLAB `atmosisa`가 가능하면 사용하고, 그렇지 않으면 fallback table / lower atmosphere model을 쓴다.

### `Reentry_Propagator.m`

120 km 이후 atmospheric entry를 전파한다. Drag, lift, heat flux, dynamic pressure, g-load, LOS geometry를 모두 기록한다.

### `J2PolarHohmann.py`

Python-side propagation model과 burn model helper가 들어 있다.

### `J2PolarHohmannShooting.py`

Phase 1 parameter optimization과 MATLAB JSON export를 담당한다.

## 11. 권장 workflow

정밀 run을 하려면 다음 순서를 권장한다.

1. Python shooting으로 Phase 1 값을 찾는다.
2. `configs/latest_python_solution.json`을 저장한다.
3. MATLAB에서 `RENDEZVOUS_CONFIG_JSON`으로 JSON을 불러온다.
4. `RENDEZVOUS_PHASE3_MODE`로 Phase 3 mode를 선택한다.
5. `RENDEZVOUS_REENTRY_SHAPE`로 shape를 선택한다.
6. MATLAB run 결과의 Phase 1 miss, Phase 3 FPA, Phase 4 heat/LOS diagnostic을 확인한다.
7. 필요한 경우 Python optimization 조건이나 JSON re-entry option을 조정한다.

## 12. 현재 limitation

- Re-entry `Cd`는 constant assumption이다.
- L/D curve는 PPT image를 바탕으로 한 approximate digitization이다.
- Full ECEF dynamics / Earth orientation model은 없다.
- Atmosphere는 ISA76-style이며, space weather / density variability는 없다.
- Plasma blackout은 모델링하지 않았다.
- Phase 2 guidance는 waypoint-impulsive이며 navigation filter나 closed-loop controller가 아니다.
- `FINITE_BURN`은 fixed-direction burn이며, thrust vector steering이나 attitude control coupling은 단순화되어 있다.

---

# 결론 및 제언

현재 코드는 mission-level rendezvous-to-entry simulator로서 구조가 꽤 명확해졌다. Python은 Phase 1 parameter search와 JSON export를 담당하고, MATLAB은 mission propagation과 visualization을 담당한다. Active MATLAB path에서는 방치된 Lambert option과 legacy Phase 3 block을 제거했고, finite burn / atmospheric drag / re-entry / LOS diagnostic이 같은 workflow 안에서 동작한다.

다음 단계로 가장 추천하는 개선은 세 가지이다.

1. **Re-entry aerodynamic database 강화**
   현재 PPT 기반 L/D와 constant Cd를 사용하고 있으므로, AoA-Mach-altitude별 `Cd`, `Cl`, `Cm` table이 있으면 trajectory와 heat/g-load prediction이 훨씬 설득력 있어진다.

2. **Communication blackout model 추가**
   지금 LOS는 geometry-only이다. Blackout-free communication을 목표로 한다면 plasma frequency, electron density correlation, RF frequency, link margin을 추가해야 한다.

3. **Regression test / batch runner 추가**
   기능이 많아졌으므로, 대표 scenario를 자동으로 돌려 Phase 1 miss, Phase 3 FPA, Phase 4 max heat flux, LOS maintained 여부를 저장하는 batch regression script가 필요하다.

---

# 참고문헌

1. COESA, NOAA, NASA, U.S. Air Force, [U.S. Standard Atmosphere, 1976](https://ntrs.nasa.gov/citations/19770009539), NASA Technical Reports Server.
2. COESA, NOAA, NASA, U.S. Air Force, [U.S. Standard Atmosphere, 1976 PDF](https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf), NOAA archive.
3. K. Sutton and R. A. Graves Jr., [A General Stagnation-Point Convective-Heating Equation for Arbitrary Gas Mixtures](https://ntrs.nasa.gov/api/citations/19720003329/downloads/19720003329.pdf), NASA TR R-376, 1971.
4. W. H. Clohessy and R. S. Wiltshire, [Terminal Guidance System for Satellite Rendezvous](https://arc.aiaa.org/doi/10.2514/8.8704), Journal of the Aerospace Sciences, 1960.
5. D. A. Vallado, *Fundamentals of Astrodynamics and Applications*, Microcosm Press, 4th ed.
6. R. R. Bate, D. D. Mueller, and J. E. White, *Fundamentals of Astrodynamics*, Dover Publications.
7. 사용자 제공 자료, `4가지 최적 형상 및 형상에 대한 LD 그래프.pptx`, re-entry vehicle geometry and L/D trend source.
