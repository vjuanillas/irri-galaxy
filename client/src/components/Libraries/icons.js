import { library } from "@fortawesome/fontawesome-svg-core";
import { faCheckSquare, faFile, faFolder, faSave, faSquare } from "@fortawesome/free-regular-svg-icons";
import {
    faAngleDoubleLeft,
    faBan,
    faBook,
    faDownload,
    faGlobe,
    faHome,
    faKey,
    faMinusSquare,
    faPencilAlt,
    faPlus,
    faShieldAlt,
    faSpinner,
    faTimes,
    faTrash,
    faUnlock,
    faUsers,
} from "@fortawesome/free-solid-svg-icons";

const tableIcons = [
    faFile,
    faFolder,
    faSpinner,
    faShieldAlt,
    faKey,
    faGlobe,
    faSave,
    faTimes,
    faBan,
    faUnlock,
    faPencilAlt,
    faUsers,
    faCheckSquare,
    faSquare,
    faMinusSquare,
];

const manageIcons = [faAngleDoubleLeft, faSave, faFile];
const topBarIcons = [faHome, faPlus, faTrash, faDownload, faBook];
const librariesIcons = [faGlobe, faPencilAlt, faSave, faTimes, faTrash, faUsers, faHome, faUnlock];

export function initFolderTableIcons() {
    tableIcons.forEach((icon) => library.add(icon));
}

export function initPermissionsIcons() {
    manageIcons.forEach((icon) => library.add(icon));
}
export function initLibrariesIcons() {
    librariesIcons.forEach((icon) => library.add(icon));
}
export function initTopBarIcons() {
    topBarIcons.forEach((icon) => library.add(icon));
}
