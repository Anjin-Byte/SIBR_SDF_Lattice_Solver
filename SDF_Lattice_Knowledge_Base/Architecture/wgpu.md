---
title: wgpu
tags: [reference, dependency, gpu]
type: external-dependency
---

# wgpu

A native Rust GPU API based on the WebGPU specification, running on Vulkan, Metal, D3D12, and OpenGL on native platforms.

## Mental model for this project

> A **portable GPU layer** for a native systems app — not a web framework, not an engine.

The "WebGPU" name can mislead. wgpu is a native Rust library that happens to expose a WebGPU-shaped API. On macOS it talks to Metal, on Windows D3D12, on Linux Vulkan. The browser is just one of its many targets.

## Why we use it

- Cross-platform GPU access without writing per-backend code
- Modern, validation-rich API (helpful during development)
- Active Rust ecosystem
- WGSL shaders compile once and run everywhere wgpu runs

## Where it lives in the workspace

**Only** in [[gpu crate]]. See [[GPU Boundary]] for the rationale and [[Shared Dependencies]] for how the version is pinned.

## Lifecycle inside `gpu`

Standard flow:
1. `Instance` — backend selection
2. `Adapter` — physical device pick
3. `Device` + `Queue` — logical device, command submission
4. `ShaderModule` from WGSL files under `crates/gpu/shaders/`
5. `BindGroupLayout` + `PipelineLayout` + `ComputePipeline` / `RenderPipeline`
6. `Buffer` / `Texture` allocation, upload via `Queue::write_buffer` or staging
7. `CommandEncoder` → `Queue::submit`
8. Optional readback via mapped staging buffer

For headless compute (the default for the SDF solver), no `Surface` is needed — see [[Headless First]].

For preview, a surface is created but lives in [[preview crate]], which asks `gpu` for an appropriate device/format pair.

## Things to remember

- Validation messages from wgpu are extremely helpful — wire them up via `tracing`.
- Use debug labels (`label: Some("...")`) on every resource. The cost is nothing; the payoff during graphics debugging is enormous.
- Pin the version in `[workspace.dependencies]` and bump deliberately — wgpu has breaking changes between minor versions historically.

## Related

- [[gpu crate]]
- [[GPU Boundary]]
- [[winit]] — companion library for windowing, intentionally separate
