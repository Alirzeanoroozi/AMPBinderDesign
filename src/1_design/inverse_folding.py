import argparse

import torch
from pytorch_lightning import Trainer
torch.set_grad_enabled(False)

from boltzgen.task.predict.data_from_generated import FromGeneratedDataModule, DataConfig
from boltzgen.data.tokenize.tokenizer import Tokenizer
from boltzgen.data.feature.featurizer import Featurizer
from boltzgen.task.predict.writer import DesignWriter
from boltzgen.model.models.boltz import Boltz

MOL_DIR = "/home/anoroozi25/.cache/huggingface/hub/datasets--boltzgen--inference-data/snapshots/c3d36fd276e9caf098c75d4113c6d5eb320b1a4c/mols.zip"
CHECKPOINT = "/home/anoroozi25/.cache/huggingface/hub/models--boltzgen--boltzgen-1/snapshots/c1be29e1f82ffcc72264f64b993c43fb4e0d17f0/boltzgen1_ifold.ckpt"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--design_dir", type=str, default="outputs/test_design_rbx1_mine/intermediate_designs")
    parser.add_argument("--output_dir", type=str, default="outputs/test_design_rbx1_mine/inverse_folding/")
    parser.add_argument("--moldir", type=str, default=MOL_DIR)
    parser.add_argument("--checkpoint", type=str, default=CHECKPOINT)
    parser.add_argument("--multiplicity", type=int, default=10)
    parser.add_argument("--samples_per_target", type=int, default=1000000000)
    parser.add_argument("--num_workers", type=int, default=1)
    parser.add_argument("--skip_existing", action="store_true")
    args = parser.parse_args()

    data_config = DataConfig(
        moldir=args.moldir,
        multiplicity=args.multiplicity,
        tokenizer=Tokenizer(atomize_modified_residues=False,),
        featurizer=Featurizer(), # TODO: Add featurizer from config
        suffix=".cif",
        suffix_metadata=".npz",
        suffix_native="_native.cif",
        samples_per_target=args.samples_per_target,
        design=True,
        backbone_only=True,
        atom14=False,
        max_seqs=1,
        inverse_fold=True,
        batch_size=1,
        num_workers=args.num_workers,
        pin_memory=True,
        num_targets=1000000000,
        design_mask_override=None
    )

    data = FromGeneratedDataModule(
        data_config,
        fail_if_no_designs=True,
        design_dir=args.design_dir,
        output_dir=args.output_dir,
        skip_existing=args.skip_existing,
        skip_existing_kind="inverse_fold",
    )

    writer = DesignWriter(
        output_dir=args.output_dir,
        res_atoms_only=False,
        atom14=data_config.atom14,
        atom37=False,
        backbone_only=False,
        inverse_fold=True,
        write_native=False,
    )

    predict_args = {
        "recycling_steps": 3,
        "sampling_steps": 200,
        "diffusion_samples": 10,
    }

    override = {
        "masker_args": {
            "mask": True,
            "mask_backbone": False,
        },
        "validators": None,
        "diffusion_process_args": {
            "sigma_min": 0.0004,
            "sigma_max": 160.0,
            "sigma_data": 16.0,
            "rho": 7,
            "P_mean": -1.2,
            "P_std": 1.5,
            "gamma_0": 0.8,
            "gamma_min": 1.0,
            "noise_scale": 1.0,
            "step_scale": 1.0,
            "mse_rotational_alignment": True,
            "coordinate_augmentation": True,
            "alignment_reverse_diff": True,
            "synchronize_sigmas": False,
        },
        "inverse_fold_args": {
            "atom_s": 128,
            "atom_z": 16,
            "token_s": 384,
            "token_z": 128,
            "node_dim": 128,
            "pair_dim": 128,
            "hidden_dim": 128,
            "dropout": 0.1,
            "softmax_dropout": 0.2,
            "num_encoder_layers": 6,
            "transformation_scale_factor": 1.0,
            "inverse_fold_noise": 0.2,
            "topk": 30,
            "num_heads": 4,
            "num_decoder_layers": 3,
            "autoregressive": True,
            "enable_input_embedder": True,
            "inverse_fold_restriction": ["CYS"],
        },
        "use_kernels": True,
    }

    # Load model
    model_module = Boltz.load_from_checkpoint(
        args.checkpoint,
        strict=True,
        use_ema=False,
        checkpoint_diffusion_conditioning=False,
        map_location="cpu",
        predict_args=predict_args,
        weights_only=False,
        **override,
    )
    model_module.eval()

    lightning_trainer = Trainer(
        default_root_dir=args.output_dir,
        strategy="auto",
        callbacks=writer,
        logger=False,
        accelerator="gpu",
        devices=1,
        precision="bf16-mixed",
    )

    # Run training
    lightning_trainer.predict(model_module, datamodule=data)
