# backend/engine.py

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

class BlackPhysicsEngine:
    def calculate(self, mission_id: str, smiles: str, user_params: dict = None) -> dict:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"pass": False, "error": "Invalid SMILES", "val": 0}

        # 基礎データ抽出
        mw = Descriptors.MolWt(mol)
        
        # 0除算ガード：分子量が極端に小さい（または0）場合は計算を中断し、事実を返す
        if mw < 0.1:
            return {
                "pass": False, 
                "error": f"物理的破綻: 分子量({mw})が小さすぎます。実在可能な構造ではありません。", 
                "val": 0,
                "smiles": smiles
            }

        mr = Crippen.MolMR(mol)
        psa = Descriptors.TPSA(mol)
        rings = Descriptors.NumAromaticRings(mol)
        num_atoms = mol.GetNumAtoms()
        rot_bonds = Descriptors.NumRotatableBonds(mol)

        res = {"smiles": smiles}
        success = False

        if mission_id == "m1":  # Stealth
            val = (mr / mw) * 100
            label, unit, target = "NIR反射率", "%", float(user_params.get("nir_reflectance", 1.5))
            success = val < target
        elif mission_id == "m2":  # Void
            val = (psa / (mw * rings + 1)) * 1000000 
            label, unit, target = "Visible_Reflectance", "ppm", float(user_params.get("visible_reflectance", 350.0))
            success = val < target # ターゲット(350)より低い数値なら成功
        elif mission_id == "m3":  # Thermal
            val = (psa / mw) * (1 + 0.05 * rot_bonds)
            label, unit, target = "放射率 ε", "mε", float(user_params.get("emissivity", 970))
            success = val < target
        elif mission_id == "m4":  # Meta
            val = 1.0 - (rings / (num_atoms + 1)) * 2.0
            label, unit, target = "負の屈折率 n", "n", float(user_params.get("refractive_index", -1.0))
            success = val < target
        else:
            return {"pass": False, "error": "Unknown Mission"}

        res.update({"val": round(val, 4), "label": label, "unit": unit, "target": target, "pass": success})
        return res
    
    # 既存のcalculateメソッドに追加
    def calculate_physics(self, mission_id: str, smiles: str):
        # すでに作成されているcalculateメソッドを呼び出すだけ
        return self.calculate(mission_id, smiles, user_params={})