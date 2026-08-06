document.addEventListener("DOMContentLoaded", () => {
    const statusEl = document.getElementById("status");
    const versionEl = document.getElementById("version");

    // createChapro is exported by the Emscripten MODULARIZE option
    if (typeof createChapro === 'function') {
        createChapro().then((Module) => {
            statusEl.textContent = "WebAssembly module loaded successfully!";
            statusEl.style.backgroundColor = "#e8f5e9";
            statusEl.style.color = "#2e7d32";

            // Bind the cha_version C function
            const cha_version = Module.cwrap('cha_version', 'string', []);
            
            // Call the function and update the DOM
            versionEl.textContent = cha_version();
        }).catch((err) => {
            statusEl.textContent = "Error loading WebAssembly module.";
            statusEl.style.backgroundColor = "#ffebee";
            statusEl.style.color = "#c62828";
            console.error(err);
        });
    } else {
        statusEl.textContent = "chapro.js not found. Make sure to build the Wasm module first.";
        statusEl.style.backgroundColor = "#fff3e0";
        statusEl.style.color = "#ef6c00";
    }
});
