<template>
    <div ref="target" />
</template>

<script>
import GridView from "mvc/grid/grid-view";
import { getAppRoot } from "onload";

export default {
    props: {
        urlBase: {
            type: String,
            default: null,
        },
        userFilter: {
            type: String,
            default: null,
        },
    },
    mounted() {
        new GridView({
            active_tab: "shared",
            url_base: `${getAppRoot()}${this.urlBase}`,
            url_data: {
                "f-username": this.userFilter || "",
            },
        }).$el.appendTo(this.$refs.target);
    },
};
</script>
