---
title: Black Core Logic Gate
emoji: ⚙️
colorFrom: gray
colorTo: gray
sdk: docker
pinned: false
---

<p align="left">
  <img src="BCLG1.png" height="150" alt="App Screenshot 1" />
  <img src="BCLG2.png" height="150" alt="App Screenshot 2" />
</p>

# Black-Core: Molecular Exploration Logic Gate

**[Available in English & Japanese]** 

A logic-driven multi-agent framework for autonomous molecular discovery and entropy analysis with LangGraph and RDkit. Available on Hugging Face: [Black Core Logic Gate](https://huggingface.co/spaces/ht55-hf/Black_Core_Logic_Gate) (Runs in a Docker container)

---

## Why This Exists

I've always been drawn to black — not as an absence of color, but as a physical phenomenon with measurable, engineerable properties. Vantablack absorbs 99.965% of light. The F-35's surface coating makes an aircraft invisible to radar. There's something deeply compelling about a color that can do that.

The deeper motivation was about LLMs. I wanted to see where they break under hard scientific constraints, where hallucinations aren't just wrong but specifically, chemically wrong in ways that matter. The answer was to pair LLM creativity with RDKit as a deterministic gatekeeper. Chemical validity is binary: a molecule either passes or it doesn't, and that hard constraint turned out to be what makes the LLM's creativity useful rather than dangerous.
Serendipity mode exists because sometimes a structurally questionable output points in a direction no one thought to look. For researchers, that's not a bug.

---

## System Architecture: The Logic-Driven Pipeline

This system operates through a series of specialized nodes and routers to balance creative molecular mutation with rigorous chemical validation.

---
> **License:** Polyform Noncommercial — research and personal use only.
> Commercial use is not permitted without explicit permission.
---

### 1. Routing & Orchestration

**Design Router:** Manages the flow between different exploration protocols (Standard, Taboo, Serendipity). It ensures the LLM's creativity is channeled into the correct mission constraints.

**Physics Router:** Bridges the non-deterministic output of LLMs with the deterministic world of cheminformatics. It handles the handover of SMILES strings to the RDKit engine for objective property calculation.

### 2. The Nodes

**Design Node (Vermouth):** Generates molecular candidates based on the mission's entropy targets.

**Mutate Node (The "Accident" Engine):**
Crucial Design Choice: Unlike traditional deterministic algorithms using RDKit for structural mutation, this node intentionally utilizes LLM-driven stochastic　mutation. By doing so, it introduces "controlled accidents" — structural transitions that transcend human-programmed heuristics — while remaining bounded by subsequent physical feasibility checks.

**Physics Node (RDKit Core):** The "Gatekeeper." It performs Kekulize validation, 3D conformation generation (SDF), and calculates molecular descriptors (MW, PSA, Rings). This node effectively filters out "hallucinated" molecules.

**Analysis Node (Sherry):** Performs a post-hoc logical evaluation of the physics data to decide the next iteration or archive the result.

**Robust Anti-Hallucination Framework**
Each LLM-powered node is fortified with **multiple layers of anti-hallucination protocols**. Every line of code and every prompt directive has been meticulously engineered to suppress the stochastic 'fantasies' of the model. By combining rigorous state-machine constraints, prompt-based role-consistency enforcement, and downstream physical validation (RDKit), the system ensures that the agents' outputs remain grounded in chemical reality and logical coherence.

#### [ Deep Dive: The "Logic Gate" Architecture ]

**1. Design & Mutate Node (Strategic Variation)**

This node functions as the creative engine, but its output is governed by mission-specific constraints.

- **LLM-Driven Stochastic Mutation:** Instead of relying on predefined fragment libraries, the node uses LLMs to perform "soft-logic" mutations. This allows for structural leaps that maintain a 2D pharmacophore's essence while introducing novel topology.

- **Recursive Feedback Loop:** The node consumes the previous iteration's distilled_lesson (from the Analysis Node) to avoid repetitive structural failures, effectively creating a "short-term memory" for chemical space exploration.

**2. Physics Node (Deterministic Validation)**

The Physics Node acts as the rigid "Gatekeeper." It translates linguistic molecular descriptions into hard scientific data.

- **Kekulization & Sanity Check:** Every SMILES string undergoes a rigorous sanitization process using RDKit. If the LLM generates a chemically "illegal" structure, the node catches the error (e.g., KekulizeException) and converts it into a System Error report.

- **Structural Quantization:** It computes multi-dimensional descriptors including:

    - **Molecular Weight (MW) Delta:** Monitoring mass changes relative to the scaffold.

    - **Polar Surface Area (PSA) & Ring Count:** Ensuring the "disguise" doesn't compromise drug-likeness.

    - **2D Conformation Scoring:** Using the ETKDG method to verify if the 2D graph can exist in a low-energy 3D state.

**3. Router Logic (Asynchronous Orchestration)**

The Routers are the "traffic controllers" of the system, ensuring data integrity across the state machine.

- **Design Router:** Determines the exploration depth. It dynamically selects the agent's "tone" and "constraint set" (Standard/Taboo/Serendipity) before the Design Node is invoked.

- **Physics Router:** A critical translation layer. It extracts the raw SMILES from the agent's dialogue, cleanses markdown artifacts, and prepares the payload for the RDKit compute engine. This prevents "parsing hallucinations" from crashing the downstream analysis.

**4. Analysis Node (Logical Synthesis)**

The final gate where data meets strategy.

- **Entropy Alignment Check:** Compares the physical mutation_score against the mission's target_entropy.

- **Strategic Decision Making: Human-in-the-Loop (HITL) Integration** Based on the results, it triggers one of two logical states:

    - "NEXT ITERATION": If the target is not met, it generates a distilled_lesson for the Design Node.

    - "SAVE RESULT": If the mission goals are satisfied, it terminates the loop and locks the data.

While the engine automates molecular generation and physical validation, the final strategic trajectory is governed by a Human-in-the-Loop interface. The system presents refined candidates and scientific insights to the user, who acts as the final decision-maker, choosing to either commit to the "NEXT ITERATION" for further refinement or "SAVE RESULT" once the mission objectives are met. This ensures that the AI's exploration remains aligned with human expertise and strategic intent.

---

## Technical Architecture

- **Graph Engine (LangGraph):** Manages the iterative feedback loop between the Molecular Designer (Vermouth), the Physicist (RDKit/Engine), and the Analyst/Chemist (Sherry). Plus the Commander (Gin). *(The agents' names were inspired by the Japanese comic/anime Detective Conan's The Black Organization.)*

- **Physical Guardrail:** Every proposed structure is subjected to 3D conformation embedding and energy minimization using RDKit to ensure structural realism before property calculation.

- **Entropy Control:** A user-defined "Entropy Slider" (0 - 100) dictates the intensity of structural transformation, allowing researchers to choose between incremental improvement and disruptive innovation.

---

## Technical Implementation (Cheminformatics Integration)

The engine ensures high-fidelity results by integrating industrial-grade tools directly into the agentic loop:

**Molecular Validation:** Every candidate is passed through RDKit to ensure structural validity.

**2D & 3D Conformation Engine:** Real-time generation of 2D coordinates from SMILES to provide spatial insights. Interactive 3D base-molecular rendering using 3Dmol.js.

**Entropy Scoring:** A custom logic that evaluates "structural disguise" by measuring the delta between candidate properties and target entropy thresholds.

---

## The Missions

Each mission utilizes the Logic Gate to navigate specific chemical constraints and physical targets, using specialized base scaffolds.

**Mission 1: Stealth Black (Radar/IR Invisibility)**

- **Objective:** To shatter the existing "1.5% reflectance" barrier of current stealth technology (e.g., F-35 aircraft) and propose the ultimate "heatless black" that neutralizes thermal detection.

- **Base Scaffolds:** Azo-Linker and Phthalocyanine.

- **Focus:** Balancing high electromagnetic absorption with low thermal emissivity.

**Mission 2: Void Black (The Light Graveyard)**

- **Objective:** To surpass the world-record "350ppm" light absorption of Vantablack, proposing a "Grave of Light" where photon escape is physically minimized.

- **Base Scaffolds:** Perylene-Core and Coronene.

- **Focus:** Maximizing structural light-trapping through extreme aromatic density.

**Mission 3: Thermal Black (Radiative Stability)**

- **Objective:** To shatter the existing 970 mε efficiency record (Stanford). "Pry open" the Atmospheric Window (8–13 μm) to eject thermal energy into the abyss of deep space.

- **Base Scaffolds:** Porphyrin-Ring and Perylene-Diimide.

- **Focus:** High-emissivity design for passive radiative cooling under extreme heat.

**Mission 4: Meta Black (The Information Abyss)**

- **Objective:** To engineer a negative refractive index (n) that warps space itself. Achieve the ultimate "Transparency of Existence" where light is not absorbed, but bypassed.

- **Base Scaffolds:** Super-Benzene and Graphene-Fragment.

- **Focus:** Non-classical photonic interaction through high-entropy topology.

---

## Mission Protocols: Core Philosophy

This engine moves away from unconstrained SMILES generation, which often leads to "hallucinated" or non-synthesizable structures. Instead, it anchors the discovery process on verified chemical scaffolds.

**Standard:** High-fidelity refinement. Reliable, stable molecular optimization. Focused on localized substitution to optimize properties without altering the core skeleton.

**Taboo:** Radical modification. High-risk exploration that defies conventional stability constraints. Permits heteroatom substitution and skeletal restructuring to break through known physical limits.

**Serendipity:** Emergent discovery with LLM. A high-entropy protocol designed to trigger "productive accidents" through the Black-Core mutation logic.

---

## Empirical Model Selection & Validation

The choice of the Gemini 2.0+ series (including Flash/Flash-Lite) was not arbitrary. Extensive comparative testing across multiple major LLM families revealed that the Gemini 2.0 architecture exhibits a superior 'Scientific Integrity' profile. In my tests, it consistently outperformed other models in following strict SMILES syntax and maintaining logical consistency under chemical constraints. This system leverages the precise instruction-following capabilities of the 2.0 series to ensure that the generative process remains grounded in empirical truth.

---

## Engineered for High-Efficiency Models

While most agentic frameworks rely on heavy, expensive LLMs, this system is uniquely engineered to thrive on lightweight, resource-constrained models (e.g., Gemini 2.5 Flash Lite). To overcome the inherent limitations of smaller models, the architecture employs:

**Granular Logic Distribution:** Complex reasoning is broken down into micro-tasks across specialized nodes, preventing model cognitive overload.

**Hard-Coded Guardrails:** Rigorous state-machine constraints and RDKit-driven validation act as 'logical exoskeletons' to support the LLM's decision-making.

**High-Density Prompt Engineering:** Every instruction is meticulously optimized for token efficiency and instruction-following precision.

This project proves that with superior architectural design, professional-grade scientific discovery is possible even within the constraints of free-tier, lightweight AI ecosystems - making it an ideal framework for individual AI developers like myself :)

---

## Tech Stack

### Backend
- Orchestration: LangGraph
- LLM: Gemini 2.0/2.5 Flash Lite
- Informatics: RDKit (Descriptors & 3D Embedding)
- API: FastAPI (Server-Sent Events)

### Frontend
- UI: React + Tailwind CSS
- Animation: Framer Motion
- 3D Rendering: 3Dmol.js

---

## Project Structure

```text
├── backend/
│   ├── main.py          # SSE streaming & API endpoints
│   ├── engine.py        # RDKit integration (Physics & SDF gen)
│   ├── graph_engine.py  # LangGraph state & agent logic
│   └── .env             # API Keys
└── frontend/
    └── src/
        ├── MissionWindow.tsx  # Primary UI & 3Dmol.js logic
        └── components/MissionHub.tsx 
```

---

## ⚖️ License
This project is licensed under the **PolyForm Noncommercial 1.0.0 License**. 
Individual, educational, and research use is highly encouraged. 
For commercial inquiries, please contact the developer.

See LICENSE for full terms.

---

## 🇯🇵 Japanese　
本プロジェクトは、複数のAIエージェントによる自律的な分子探索と、堅実なロジック・物理シミュレーションを統合した研究開発プラットフォームです。LangGraphを用いたマルチエージェント・オーケストレーションにより、特定の物理特性（近赤外反射率、負の屈折率等）を持つ黒色分子構造を動的に生成・検証します。ユーザーによるパラメータ指定（scaffold, standard/taboo philosophy, entropy level）と、AIによる偶発的突然変異（serendipity philosophy）のハイブリッドな設計フローで、LangGraphの同じグラフ構造を使いつつ、Stateにどのモードで動くかのフラグを持たせ異なる挙動をさせ、ユーザーのターゲット数値、もしくは人間の想像を超えたアクシデントを人為的にかつランダムに引き起こし、ベンタブラックやF-35戦闘機用の塗装のような既存の限界を超える黒の発見を目指すこと、の両方の作動ができるようにしました。

新しい究極の黒色分子の探究をゴールにしているため、PythonとRDkitでガチガチに縛りつけるよりも、AIの持つ人間の想像を超えるような未知の発想の可能性を引き出すために、３つのLLM介入ノードを使用しつつも(Mutate nodeにも敢えてLLMを採用した理由の１つでもあります)、それぞれのLLMの自由度を的確にコントロールし、各所に幾重にも張り巡らせたハルシネーション対策を施すことで、LLMの暴走とハルシネーションを極限まで抑えた設計を構築しました。パラメータの種類やBase ScaffoldのSMILESなどを変えることで、さらに別の角度からの新分子の発見も期待できます。黒ずくめの組織を意識したエージェント達が、４つの異なる性質の”究極の黒”を作るミッションの下に、AIと人間のInteractionの方法や深さを変えながら作ってあります。

特に本システムでは、制約の厳しい軽量モデル（Gemini 2.0/2.5 Flash Lite）の無料枠をベースに作ったので、軽量モデルの可能性を極限まで引き出したルーティング・プロンプトになっています。コードの１つ１つ、プロンプトの一言一言が、気の遠くなるような試行錯誤を経て作られています。なお、無料LLMモデルのセレクションにあたり、複数のメジャーなLLMで基礎システムを回してテストした結果、一番賢く、科学的なハルシネーション(例えば、科学者であるはずのシェリーが”科学的な事実”と称して巧妙に嘘を言うなど)が一番低かったGemini2.0以上のFlashモデルに決めたという経緯があります。

本システムにはユーザー自身のGoogle API keyが必要です。
１つのミッションにつき、最低でも３回以上叩き、出力によってはToken消費量も早いことが予想されます。無料枠をご利用の方は、ご自身のAPI制限に気を付けつつも、もしリミットに達した場合には翌日まで待って再度お試しください。

- **Standard/Taboo Mode**: ユーザーが指定したエントロピーを目標にして、AIエージェント達が科学的に妥当で、実用性の高い堅実な最適解を導き出す。
- **Serendipity Mode**: LLMお得意の狂気モード。ベルモット達に人間では思いつかない、飛躍した（時に狂気的な）新素材をぶち上げてもらう。もちろん物理ノードとシェリーのアナリシスノードでガッツリと科学的検証を行うため、時には解析不能エラーも有り。
