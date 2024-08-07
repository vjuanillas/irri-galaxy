import { BIO_COMPUTE_OBJ_EXPORT_PLUGIN } from "./BioComputeObject/BioComputeObjectExportPlugin";
import { DEFAULT_FILE_EXPORT_PLUGIN } from "./DefaultFileExportPlugin";
import { RO_CRATE_EXPORT_PLUGIN } from "./ROCrateExportPlugin";

const AVAILABLE_INVOCATION_EXPORT_PLUGINS = [
    RO_CRATE_EXPORT_PLUGIN,
    BIO_COMPUTE_OBJ_EXPORT_PLUGIN,
    DEFAULT_FILE_EXPORT_PLUGIN,
];

export { AVAILABLE_INVOCATION_EXPORT_PLUGINS };
